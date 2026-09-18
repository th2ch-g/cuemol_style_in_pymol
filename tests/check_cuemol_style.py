"""Run real PyMOL checks in a fresh process, outside the suite's PyMOL mocks.

Use --gui for OpenGL, Qt mouse events, and image generation. Output belongs in
an ignored directory. No network access or user startup files are required.
"""

import argparse
import json
import sys
from contextlib import contextmanager
from io import BytesIO
from pathlib import Path
from time import monotonic, perf_counter, sleep
from unittest.mock import patch

import numpy as np


def snapshot(cmd, selection="all"):
    rows = []
    cmd.iterate(
        selection, "rows.append((model,index,reps,color,ss))", space={"rows": rows}
    )
    return rows, cmd.get_coords(selection).copy(), cmd.get_setting_tuple("bg_rgb")


def same_snapshot(cmd, before, selection="all"):
    after = snapshot(cmd, selection)
    assert before[0] == after[0]
    np.testing.assert_array_equal(before[1], after[1])
    assert before[2] == after[2]


def pixels(data):
    from PIL import Image

    return np.asarray(Image.open(BytesIO(data)).convert("RGB"))


def native_ray(cmd, width=240, height=180):
    cmd.ray(width, height, quiet=1)
    return cmd.png(None, prior=1, quiet=1)


def compare_native_ray(cmd, manager):
    """Compare retained ray proxies with equivalent regular triangle CGO."""
    from cuemol_style_in_pymol.export import settings

    # Isolate geometry and colors from native shadow/outline sampling noise.
    with settings(cmd, {"ray_trace_mode": 0, "ray_shadows": 0, "antialias": 0}):
        return _compare_native_ray(cmd, manager)


def _compare_native_ray(cmd, manager):
    from cuemol_style_in_pymol.export import cgo_mesh

    actual = pixels(native_ray(cmd))
    assert len(np.unique(actual.reshape(-1, 3), axis=0)) > 30
    name = "_ray_reference"
    values = []
    for drawing in manager.active_drawings():
        for piece in drawing.pieces:
            values.extend(cgo_mesh(piece, drawing.profile.material))
    cmd.load_cgo(values, name, zoom=0)
    enabled = set(cmd.get_names("objects", enabled_only=1))
    disabled = [
        n for e in manager.entries.values() for n in e.generated if n in enabled
    ]
    for obj in disabled:
        cmd.disable(obj)
    try:
        expected = pixels(cmd.png(None, 240, 180, ray=1, quiet=1))
        np.testing.assert_allclose(actual, expected, atol=1)
    finally:
        cmd.delete(name)
        for obj in disabled:
            cmd.enable(obj)
    return actual


def peptide(cmd, structure=None):
    if structure:
        cmd.load(str(structure), "peptide")
    else:
        cmd.fab("ACDEFGHIKLMN", name="peptide", ss=1)
        cmd.alter("peptide and resi 8-10", 'ss="S"')
    cmd.remove("hydro")
    cmd.show_as("cartoon", "peptide")
    cmd.color("salmon", "peptide")
    cmd.bg_color("grey30")
    cmd.orient("peptide")
    cmd.zoom("peptide", 3)
    cmd.set("orthoscopic", 1)


def registration_checks(cmd, report):
    from cuemol_style_in_pymol import __init_plugin__, cuemol_style
    from cuemol_style_in_pymol.controller import manager_for

    baseline = snapshot(cmd)
    __init_plugin__()
    assert cmd.keyword["cuemol_style"][0] is cuemol_style
    same_snapshot(cmd, baseline)
    cmd.do("cuemol_style toon1, quality=low, quiet=1")
    assert "cuemol" in manager_for(cmd).entries
    cmd.do("cuemol_style reset")
    assert not manager_for(cmd).entries
    same_snapshot(cmd, baseline)
    report["command_registration"] = "standalone import and PyMOL command dispatch"


def headless_checks(cmd, output, style, report):
    from pymol import CmdException

    from cuemol_style_in_pymol import export, geometry
    from cuemol_style_in_pymol.controller import manager_for
    from cuemol_style_in_pymol.presets import PROFILES

    native_commands = (cmd.ray, cmd.png, cmd.draw, cmd.do)
    peptide(cmd)
    registration_checks(cmd, report)
    baseline = snapshot(cmd)
    profiles = []
    for profile in PROFILES:
        entry = style(profile, quality="low", quiet=1, _self=cmd)
        assert cmd.get_setting_tuple("bg_rgb") == baseline[2]
        assert all(
            np.isfinite(p.mesh.vertices).all()
            for ds in entry.drawings.values()
            for d in ds
            for p in d.pieces
        )
        compare_native_ray(cmd, manager_for(cmd))
        (output / f"native_ray_{profile}.png").write_bytes(native_ray(cmd))
        style(
            "ray",
            filename=str(output / f"headless_{profile}.png"),
            width=240,
            height=180,
            _self=cmd,
        )
        style("reset", _self=cmd)
        same_snapshot(cmd, baseline)
        profiles.append(profile)
    report["headless_profiles"] = profiles
    style("toon1", quiet=1, _self=cmd)
    manager = manager_for(cmd)
    previous = manager.entries["cuemol"]
    original_proxy = export.ray_proxy
    proxy_calls = 0

    def fail_first_proxy(drawing):
        nonlocal proxy_calls
        proxy_calls += 1
        if proxy_calls == 1:
            raise RuntimeError("injected CGO failure")
        return original_proxy(drawing)

    with patch.object(export, "ray_proxy", side_effect=fail_first_proxy):
        try:
            style("toon2", quiet=1, _self=cmd)
        except CmdException:
            pass
        else:
            raise AssertionError("A failed CGO load must raise")
    assert manager.entries["cuemol"] is previous
    compare_native_ray(cmd, manager)
    with patch.object(
        geometry, "build", side_effect=RuntimeError("injected mesh failure")
    ):
        try:
            style("toon2", quiet=1, _self=cmd)
        except CmdException:
            pass
        else:
            raise AssertionError("A failed build must raise")
    assert manager.entries["cuemol"] is previous
    enabled = set(cmd.get_names("objects", enabled_only=1))
    with patch.object(cmd, "ray", side_effect=RuntimeError("injected export failure")):
        try:
            style("ray", filename=str(output / "must_not_exist.png"), _self=cmd)
        except CmdException:
            pass
        else:
            raise AssertionError("A failed export must raise")
    assert set(cmd.get_names("objects", enabled_only=1)) == enabled
    assert not (output / "must_not_exist.png").exists()
    assert not any(n.startswith("_cuemol_ray_") for n in cmd.get_names())
    cmd.pseudoatom("cuemol_alpha_1", pos=[100, 100, 100])
    try:
        style("toon1", selection="peptide", transparency=0.4, quiet=1, _self=cmd)
    except CmdException:
        pass
    else:
        raise AssertionError("Generated names must not overwrite source objects")
    assert cmd.count_atoms("cuemol_alpha_1") == 1
    assert manager.entries["cuemol"] is previous
    assert "cuemol" in cmd.get_names("objects")
    cmd.delete("cuemol_alpha_1")
    try:
        style("matte", name="overlap", quiet=1, _self=cmd)
    except CmdException:
        pass
    else:
        raise AssertionError("Overlapping managed selections must be rejected")
    session = output / "session.pse"
    cmd.save(str(session))
    style("reset", _self=cmd)
    cmd.load(str(session))
    assert list(manager.entries) == ["cuemol"]
    compare_native_ray(cmd, manager)
    style("reset", _self=cmd)
    same_snapshot(cmd, baseline)
    style("richardson", quiet=1, _self=cmd)
    cmd.delete("cuemol_shape_1")
    manager.maintenance()
    same_snapshot(cmd, baseline)
    assert manager.entries == {}
    assert (cmd.ray, cmd.png, cmd.draw, cmd.do) == native_commands
    cmd.fragment("ala", "ligand")
    cmd.show_as("lines", "ligand")
    cmd.hide("everything", "ligand and name CB")
    mixed = snapshot(cmd)
    entry = style("toon1", representation="auto", quiet=1, _self=cmd)
    assert set(entry.drawings) == {"peptide", "ligand"}
    assert entry.drawings["peptide"][0].pieces
    atoms = [a for p in entry.drawings["ligand"][0].pieces for a in p.atoms]
    assert atoms and all(a.name != "CB" for a in atoms)
    style("reset", _self=cmd)
    same_snapshot(cmd, mixed)
    style("toon1", quiet=1, _self=cmd)
    cmd.bg_color("blue")
    manual_background = cmd.get_setting_tuple("bg_rgb")
    assert manual_background != baseline[2]
    style("refresh", _self=cmd)
    assert cmd.get_setting_tuple("bg_rgb") == manual_background
    cmd.save(str(session))
    cmd.load(str(session))
    assert cmd.get_setting_tuple("bg_rgb") == manual_background
    style("reset", _self=cmd)
    assert cmd.get_setting_tuple("bg_rgb") == manual_background
    report["headless_lifecycle"] = (
        "apply/reset, failed replacement/export, overlap, saved session, generated-object deletion, unchanged standard commands"
    )
    cmd.delete("all")
    state_checks(cmd, style)
    report["headless_states"] = (
        "nonsequential states, movie mapping, per-object overrides, all states, source enable/disable, added states, transparent state count"
    )
    report["native_ray"] = (
        "standard ray and png ray=1 match native CGO: all profiles, rotations, immediate state changes, movie mapping, mixed opacity, session reload"
    )


def state_checks(cmd, style, pump=None):
    from cuemol_style_in_pymol.controller import manager_for

    manager = manager_for(cmd)
    cmd.fragment("ala", "trajectory")
    cmd.show_as("spheres")
    base = cmd.get_coords("trajectory")
    for state in (2, 3):
        cmd.create("trajectory", "trajectory", 1, state)
        cmd.load_coords(base + [state * 4, 0, 0], "trajectory", state)
    cmd.frame(2)
    style("toon1", representation="cpk", quiet=1, _self=cmd)
    cmd.reset()
    cmd.zoom("trajectory", state=0, buffer=3)
    previous = None
    for state in (1, 3, 2):
        cmd.frame(state)
        assert cmd.get_object_state("cuemol_ray_1") == state
        current = compare_native_ray(cmd, manager)
        if previous is not None:
            assert not np.array_equal(current, previous)
        previous = current
        manager.maintenance()
        if pump:
            pump(0.12)
        assert cmd.get_object_state("cuemol_shape_1") == state
        active = list(manager.active_drawings())
        assert [d.state for d in active] == [state]
        if pump:
            assert active[0].draws > 0
    cmd.mset("1 3 2")
    cmd.frame(2)
    if pump:
        pump(0.12)
    assert cmd.get_state() == 3
    assert [d.state for d in manager.active_drawings()] == [3]
    compare_native_ray(cmd, manager)
    cmd.turn("y", 65)
    compare_native_ray(cmd, manager)
    cmd.turn("y", -65)
    cmd.set("state", 1, "trajectory")
    manager.maintenance()
    assert [d.state for d in manager.active_drawings()] == [1]
    assert cmd.get_object_state("cuemol_ray_1") == 1
    compare_native_ray(cmd, manager)
    cmd.unset("state", "trajectory")
    manager.maintenance()
    assert [d.state for d in manager.active_drawings()] == [3]
    cmd.set("all_states", 1, "trajectory")
    manager.maintenance()
    assert [d.state for d in manager.active_drawings()] == [1, 2, 3]
    compare_native_ray(cmd, manager)
    cmd.unset("all_states", "trajectory")
    cmd.disable("trajectory")
    manager.maintenance()
    assert list(manager.active_drawings()) == []
    cmd.enable("trajectory")
    manager.maintenance()
    cmd.mset("")
    cmd.create("trajectory", "trajectory", 1, 4)
    style("refresh", _self=cmd)
    assert cmd.count_states("cuemol_shape_1") == 4
    assert cmd.count_states("cuemol_ray_1") == 4
    cmd.set("sphere_transparency", 0.4, "trajectory and name N")
    style("refresh", _self=cmd)
    assert "cuemol_ray_1" in cmd.get_names("objects")
    assert "cuemol_alpha_1" in cmd.get_names("objects")
    compare_native_ray(cmd, manager)
    style("toon1", representation="cpk", transparency=0.4, quiet=1, _self=cmd)
    assert cmd.count_states("cuemol_alpha_1") == 4
    assert "cuemol_ray_1" not in cmd.get_names("objects")
    compare_native_ray(cmd, manager)
    cmd.frame(4)
    if pump:
        pump(0.12)
    assert cmd.get_object_state("cuemol_alpha_1") == 4
    style("reset", _self=cmd)
    cmd.delete("all")
    cmd.frame(1)


@contextmanager
def qt_session():
    import pymol

    pymol.invocation.options.plugins = 0
    pymol.invocation.options.deferred = []
    pymol.invocation.options.pymolrc = None
    pymol.invocation.options.no_gui = 0
    pymol.invocation.options.show_splash = 0
    from pmg_qt.pymol_qt_gui import PyMOLApplication, PyMOLQtGUI

    app = PyMOLApplication(["cuemol-style-checks"])
    window = PyMOLQtGUI()
    widget = window.pymolwidget
    cmd = widget.cmd

    def in_context(func):
        with widget:
            return func()

    # Match PyMOL's normal execapp context dispatch for this standalone harness.
    cmd._call_with_opengl_context = in_context

    def pump(seconds=0.1):
        end = monotonic() + seconds
        while monotonic() < end:
            app.processEvents()
            sleep(0.002)

    window.resize(1000, 850)
    window.show()
    pump(0.3)
    try:
        yield cmd, widget, pump
    finally:
        manager = vars(cmd._pymol).get("_cuemol_style_manager")
        if manager:
            manager.reset()
        window.hide()
        cmd.delete("all")


def gui_checks(cmd, widget, pump, output, style, report, structure):
    from OpenGL import GL as gl
    from pymol.Qt import QtCore, QtGui, QtWidgets
    from PyQt5.QtTest import QTest

    from cuemol_style_in_pymol.controller import manager_for
    from cuemol_style_in_pymol.export import view_matrix
    from cuemol_style_in_pymol.picking import hit_at
    from cuemol_style_in_pymol.presets import PROFILES

    manager = manager_for(cmd)
    with widget:
        report["opengl"] = gl.glGetString(gl.GL_VERSION).decode()
        report["renderer"] = gl.glGetString(gl.GL_RENDERER).decode()
    peptide(cmd, structure)
    registration_checks(cmd, report)
    # Interactive drawing must never re-enter the precision CPU rasterizer or
    # download full framebuffer images during rotation.
    with (
        patch(
            "cuemol_style_in_pymol.raster.buffers",
            side_effect=AssertionError("CPU contour pass"),
        ),
        patch(
            "cuemol_style_in_pymol.hatch.pencil",
            side_effect=AssertionError("CPU pencil pass"),
        ),
        patch.object(
            gl, "glReadPixels", side_effect=AssertionError("GPU image readback")
        ),
    ):
        for profile in (
            "toon1",
            "toon2",
            "outline",
            "silhouette",
            "richardson",
            "ribbon",
        ):
            entry = style(profile, quiet=1, _self=cmd)
            drawing = next(manager.active_drawings())
            count = drawing.draws
            for angle in (3, -3):
                cmd.turn("y", angle)
                with widget:
                    widget.paintGL()
                    gl.glFinish()
            assert drawing.draws >= count + 2 and not drawing.error
            assert not manager.pool.precise
            style("reset", _self=cmd)
    report["interactive_gpu"] = (
        "six styles rotate without CPU contours, pencil sampling, or image readback"
    )
    cmd.set("internal_gui", 0)
    cmd.set("internal_feedback", 0)
    baseline = snapshot(cmd)
    profiles = []
    for profile in PROFILES:
        entry = style(profile, quality="medium", quiet=1, _self=cmd)
        pump()
        style(
            "png",
            filename=str(output / f"gpu_{profile}.png"),
            width=640,
            height=480,
            _self=cmd,
        )
        assert all(
            d.draws > 0 and not d.error for ds in entry.drawings.values() for d in ds
        )
        if profile == "richardson":
            (output / "native_ray_richardson.png").write_bytes(
                native_ray(cmd, 640, 480)
            )
            shader = cmd.get_setting_int("cgo_use_shader")
            cmd.disable("cuemol_shape_1")
            try:
                for enabled_shader in (1, 0):
                    cmd.set("cgo_use_shader", enabled_shader)
                    cmd.draw(320, 240, quiet=1)
                    with_proxy = pixels(cmd.png(None, prior=1, quiet=1))
                    cmd.disable("cuemol_ray_1")
                    cmd.draw(320, 240, quiet=1)
                    without_proxy = pixels(cmd.png(None, prior=1, quiet=1))
                    np.testing.assert_array_equal(with_proxy, without_proxy)
                    cmd.enable("cuemol_ray_1")
            finally:
                cmd.set("cgo_use_shader", shader)
                cmd.enable("cuemol_shape_1")
            style(
                "png",
                filename=str(output / "gpu_high_resolution.png"),
                width=2400,
                height=1800,
                _self=cmd,
            )
            style(
                "ray",
                filename=str(output / "ray_richardson.png"),
                width=640,
                height=480,
                _self=cmd,
            )
            cmd.turn("y", 65)
            pump()
            style(
                "png",
                filename=str(output / "gpu_richardson_rotated.png"),
                width=640,
                height=480,
                _self=cmd,
            )
            cmd.turn("y", -65)
        style("reset", _self=cmd)
        same_snapshot(cmd, baseline)
        profiles.append(profile)
    report["gpu_profiles"] = profiles
    cmd.delete("all")
    cmd.pseudoatom("ligand", name="C1", resi="1", resn="LIG", elem="C", pos=[-2, 0, 0])
    cmd.pseudoatom("ligand", name="O1", resi="1", resn="LIG", elem="O", pos=[-2, 1, 0])
    cmd.pseudoatom("ligand", name="N1", resi="2", resn="LIG", elem="N", pos=[2, 0, 0])
    cmd.show_as("spheres", "ligand")
    cmd.color("salmon", "ligand")
    cmd.reset()
    cmd.zoom("ligand", 3)
    entry = style("toon1", representation="cpk", quiet=1, _self=cmd)
    pump(0.3)
    drawing = next(manager.active_drawings())
    np.testing.assert_allclose(drawing.matrices[0], view_matrix(cmd), atol=1e-5)

    def point(world):
        d = next(manager.active_drawings())
        mv, projection, viewport = d.matrices
        clip = projection @ mv @ [*world, 1]
        ndc = clip[:3] / clip[3]
        x = viewport[0] + (ndc[0] + 1) * viewport[2] / 2
        y = viewport[1] + (ndc[1] + 1) * viewport[3] / 2
        return x, y

    def click(world, shift=False):
        x, y = point(world)
        ratio = widget.devicePixelRatioF()
        pos = QtCore.QPoint(round(x / ratio), round(widget.height() - y / ratio))
        QTest.mouseClick(
            widget,
            QtCore.Qt.LeftButton,
            QtCore.Qt.ShiftModifier if shift else QtCore.Qt.NoModifier,
            pos,
        )
        pump(0.2)

    cmd.set("mouse_selection_mode", 0)
    click([-2, 0, 0])
    assert cmd.count_atoms("sele") == 1
    assert cmd.count_atoms("sele and name C1") == 1
    assert len(next(manager.active_drawings()).selected) == 1
    widget.grabFramebuffer().save(str(output / "selection_live.png"))
    cmd.set("mouse_selection_mode", 1)
    click([-2, 0, 0])
    assert cmd.count_atoms("sele") == 2
    click([2, 0, 0], shift=True)
    assert cmd.count_atoms("sele") == 3
    before = np.asarray(cmd.get_view())
    x, y = point([-2, 0, 0])
    start = QtCore.QPoint(
        round(x / widget.devicePixelRatioF()),
        round(widget.height() - y / widget.devicePixelRatioF()),
    )
    with patch.object(
        manager.events, "hit", side_effect=AssertionError("Drag picking")
    ) as pick:
        QTest.mousePress(widget, QtCore.Qt.LeftButton, pos=start)
        QtWidgets.QApplication.sendEvent(
            widget,
            QtGui.QMouseEvent(
                QtCore.QEvent.MouseMove,
                QtCore.QPointF(start + QtCore.QPoint(90, 40)),
                QtCore.Qt.NoButton,
                QtCore.Qt.LeftButton,
                QtCore.Qt.NoModifier,
            ),
        )
        assert manager.pool.interacting
        with widget:
            widget.paintGL()
            gl.glFinish()
        assert manager.pool.preview.bytes > 0
        QTest.mouseRelease(
            widget, QtCore.Qt.LeftButton, pos=start + QtCore.QPoint(90, 40)
        )
        pick.assert_not_called()
    pump(0.3)
    assert not manager.pool.interacting
    assert not np.allclose(cmd.get_view(), before)
    cmd.set_view(before.tolist())
    cmd.deselect()
    pump()
    # An opaque native sphere in front must block the custom geometry's picker.
    cmd.pseudoatom("occluder", pos=[-2, 0, 4], vdw=2.5)
    cmd.show_as("spheres", "occluder")
    cmd.color("blue", "occluder")
    pump()
    x, y = point([-2, 0, 0])
    depth = manager.events.depth(x, y)
    assert hit_at(list(manager.active_drawings()), x, y, depth) is None
    style(
        "png",
        filename=str(output / "native_occlusion.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    cmd.delete("occluder")
    style("reset", _self=cmd)
    style("toon1", representation="cpk", transparency=0.45, quiet=1, _self=cmd)
    cmd.pseudoatom("native_alpha", pos=[0, 0, 1], vdw=2.7)
    cmd.show_as("spheres", "native_alpha")
    cmd.color("cyan", "native_alpha")
    cmd.set("sphere_transparency", 0.5, "native_alpha")
    pump()
    style(
        "png",
        filename=str(output / "transparency_native.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    cmd.delete("native_alpha")
    style("reset", _self=cmd)
    # Two custom surfaces use PyMOL's native triangle transparency pass together.
    cmd.create("copy", "ligand")
    cmd.translate([1, 0, 2], "copy")
    cmd.zoom("ligand or copy", 3)
    style(
        "toon1",
        selection="ligand",
        representation="surface",
        transparency=0.4,
        name="front_view",
        quiet=1,
        _self=cmd,
    )
    style(
        "matte",
        selection="copy",
        representation="surface",
        transparency=0.55,
        name="back_view",
        color="keep",
        quiet=1,
        _self=cmd,
    )
    cmd.color("blue", "copy")
    style("refresh", name="back_view", _self=cmd)
    pump()
    assert set(manager.entries) == {"front_view", "back_view"}
    assert len(list(manager.active_drawings())) == 2
    style(
        "png",
        filename=str(output / "transparency_surfaces.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    style(
        "ray",
        filename=str(output / "transparency_surfaces_ray.png"),
        width=640,
        height=480,
        _self=cmd,
    )
    style("reset", name="all", _self=cmd)
    cmd.delete("all")
    cmd.fnab("ATGCA", name="dna")
    cmd.show_as("cartoon")
    cmd.orient()
    cmd.zoom("dna", 4)
    style("nucleic", quiet=1, _self=cmd)
    pump()
    style("png", filename=str(output / "nucleic.png"), width=640, height=480, _self=cmd)
    style("reset", _self=cmd)
    cmd.delete("all")
    report["gui_interaction"] = (
        "atom/residue/shift selection, native drag rotation, native depth occlusion, native/custom transparency, two transparent surfaces, nucleic base-pair rods"
    )
    state_checks(cmd, style, pump)
    report["gui_states"] = "native callback and CGO state switching plus movie mapping"
    cmd.fragment("ala", "peptide")
    cmd.show_as("spheres")
    baseline = snapshot(cmd)
    entry = style("toon1", quiet=1, _self=cmd)
    pump()
    keys = (
        gl.GL_CURRENT_PROGRAM,
        gl.GL_ARRAY_BUFFER_BINDING,
        gl.GL_ELEMENT_ARRAY_BUFFER_BINDING,
        gl.GL_CLIENT_ACTIVE_TEXTURE,
        gl.GL_FRAMEBUFFER_BINDING,
        gl.GL_DEPTH_WRITEMASK,
        gl.GL_DEPTH_FUNC,
        gl.GL_BLEND,
        gl.GL_CULL_FACE,
        gl.GL_FOG,
    )
    with widget:
        before = [int(gl.glGetIntegerv(key)) for key in keys]
        manager.pool.draw(entry.drawings["peptide"][0])
        assert before == [int(gl.glGetIntegerv(key)) for key in keys]
        assert gl.glGetError() == gl.GL_NO_ERROR
    style("reset", _self=cmd)
    with patch.object(
        manager.pool, "program", side_effect=RuntimeError("injected shader failure")
    ):
        style("toon1", quiet=1, _self=cmd)
        pump(0.3)
    assert not manager.entries
    same_snapshot(cmd, baseline)
    cmd.delete("all")
    report["gui_gl_state"] = (
        "GL state preserved; native view restored after shader failure"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--benchmark-only", action="store_true")
    parser.add_argument("--benchmark-style", default="richardson")
    parser.add_argument(
        "--installed",
        action="store_true",
        help="Validate the installed distribution instead of the source tree",
    )
    parser.add_argument("--structure", type=Path)
    parser.add_argument("--output", type=Path, default=Path(".cache/cuemol-checks"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if not args.installed:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    try:
        import pymol
    except ImportError:
        print("PyMOL is not installed")
        return 77
    report = {"pymol": pymol.cmd.get_version()[0], "python": sys.version.split()[0]}
    start = perf_counter()
    if args.gui:
        with qt_session() as (cmd, widget, pump):
            from cuemol_style_in_pymol import cuemol_style

            if not args.benchmark_only:
                gui_checks(
                    cmd, widget, pump, args.output, cuemol_style, report, args.structure
                )
            if args.benchmark or args.benchmark_only:
                from cuemol_benchmark import benchmark

                benchmark(
                    cmd,
                    widget,
                    pump,
                    args.output,
                    cuemol_style,
                    report,
                    args.structure,
                    args.benchmark_style,
                )
    else:
        from cuemol_ray_checks import compare_ray_paths

        from cuemol_style_in_pymol import cuemol_style

        headless_checks(pymol.cmd, args.output, cuemol_style, report)
        report["ray_sample_composition"] = compare_ray_paths(
            pymol.cmd, args.output / "ray-composition"
        )
    report["elapsed_seconds"] = perf_counter() - start
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
