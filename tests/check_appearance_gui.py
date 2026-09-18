"""Regression checks for tiled pencil output and camera-dependent native CGO."""

import argparse
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from check_cuemol_style import pixels, qt_session, same_snapshot, snapshot

from cuemol_style_in_pymol import cuemol_style
from cuemol_style_in_pymol.controller import manager_for


def check_live_pencil(cmd, widget, pump):
    """Compare GPU body pixels with the independently tested native pencil sampler."""
    from OpenGL import GL as gl
    from OpenGL.GL.EXT import framebuffer_object as fb

    from cuemol_style_in_pymol._edtsurf import pencil
    from cuemol_style_in_pymol.gpu import Pool
    from cuemol_style_in_pymol.materials import drawing_tone

    cmd.delete("all")
    cmd.frame(1)
    cmd.pseudoatom("pencil", elem="C", vdw=1, state=1)
    cmd.set_color("pencil_pigment", [0.25, 0.5, 0.75])
    cmd.color("pencil_pigment", "pencil")
    cmd.orient("pencil")
    cmd.zoom("pencil", 1)
    cmd.set("orthoscopic", 1)
    entry = cuemol_style(
        "richardson", representation="cpk", color="keep", quiet=1, _self=cmd
    )
    pump()
    drawing = entry.drawings["pencil"][0]
    pool = manager_for(cmd).pool
    with widget:
        widget.paintGL()
        if b"GL_EXT_gpu_shader4" not in gl.glGetString(gl.GL_EXTENSIONS).split():
            cuemol_style("reset", quiet=1, _self=cmd)
            return None
        width, height = pool.live.size
        framebuffer = int(gl.glGetIntegerv(gl.GL_FRAMEBUFFER_BINDING))
        read_buffer = int(gl.glGetIntegerv(gl.GL_READ_BUFFER))
        try:
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, pool.live.framebuffer)
            arrays = []
            for attachment in (gl.GL_COLOR_ATTACHMENT0, gl.GL_COLOR_ATTACHMENT1):
                gl.glReadBuffer(attachment)
                arrays.append(
                    np.asarray(
                        gl.glReadPixels(0, 0, width, height, gl.GL_RGBA, gl.GL_FLOAT)
                    ).reshape(height, width, 4)
                )
        finally:
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, framebuffer)
            gl.glReadBuffer(read_buffer)
        fallback = Pool()
        try:
            with patch.object(gl, "glGetString", return_value=b""):
                assert fallback.program("body")
        finally:
            for program in fallback.programs.values():
                gl.glDeleteProgram(program)
    actual, surface = arrays
    mask = surface[:, :, 3] < 0
    yy, xx = np.nonzero(mask)
    assert len(xx) > 1000
    fog = np.clip(
        (drawing.fog[1] + surface[mask, 3]) / (drawing.fog[1] - drawing.fog[0]), 0, 1
    )
    tones = drawing_tone(surface[mask, :3], fog=fog)
    colors = np.tile(drawing.pieces[0].mesh.colors[0], (len(xx), 1))
    expected = pencil(
        np.c_[(xx + 0.5) / 3, (height - yy - 0.5) / 3], tones, colors, False
    )
    difference = np.abs(actual[mask, :3] - expected) * 255
    error = float(difference.mean())
    assert error < 0.3, error  # RGBA8 quantization contributes about 0.25/255.
    assert np.percentile(difference, 99.9) < 1.0
    cuemol_style("reset", quiet=1, _self=cmd)
    return error


def check(output):
    from OpenGL import GL as gl
    from pymol import CmdException

    output.mkdir(parents=True, exist_ok=True)
    with qt_session() as (cmd, widget, pump):
        for setting in (
            "internal_gui",
            "internal_feedback",
            "seq_view",
            "movie_panel",
            "internal_prompt",
        ):
            cmd.set(setting, 0)
        ratio = widget.devicePixelRatioF()
        widget.setFixedSize(round(320 / ratio), round(240 / ratio))
        pump()
        with widget:
            widget.resizeGL(widget.width(), widget.height())
        cmd.fab("AAAAAA", name="protein", ss=1)
        cmd.remove("hydro")
        cmd.show_as("cartoon")
        cmd.orient()
        cmd.zoom(buffer=3)
        cmd.clip("slab", 150)
        cmd.bg_color("white")
        baseline = snapshot(cmd)
        manager = manager_for(cmd)
        # The first draw can precede the Qt maintenance timer on a large display.
        widget.setFixedSize(round(3200 / ratio), round(1800 / ratio))
        pump()
        with widget:
            widget.resizeGL(widget.width(), widget.height())
        cmd.set("orthoscopic", 0)
        with patch.object(manager, "prepare_view"):
            initial = cuemol_style("richardson", quiet=1, _self=cmd)
            with widget:
                widget.paintGL()
            drawing = initial.drawings["protein"][0]
            assert drawing.draws > 0 and not drawing.error, drawing.error
            assert tuple(drawing.matrices[2][2:]) == (3200, 1800)
            assert manager.pool.live.bytes < 64 * 1024**2
            assert drawing.fog[0] > 0
        cuemol_style("reset", quiet=1, _self=cmd)
        cmd.set("orthoscopic", 1)
        widget.setFixedSize(round(320 / ratio), round(240 / ratio))
        pump()
        with widget:
            widget.resizeGL(widget.width(), widget.height())
        entry = cuemol_style("richardson", quiet=1, _self=cmd)
        pump()
        live_errors = []
        for orthoscopic in (0, 1):
            cmd.set("orthoscopic", orthoscopic)
            images = []
            for tile in (510, 96):
                with patch.object(manager.pool.live, "tile_size", tile):
                    cmd.turn("y", 0)
                    with widget:
                        widget.paintGL()
                    image_path = output / f"live-{orthoscopic}-{tile}.png"
                    assert widget.grabFramebuffer().save(str(image_path))
                    images.append(pixels(image_path.read_bytes()).astype(float))
            live_error = float(np.abs(images[0] - images[1]).mean())
            assert live_error < 0.1, live_error
            live_errors.append(live_error)
        for name, tile in (("single", 510), ("tiled", 96)):
            with patch.object(manager.pool.hatch, "tile_size", tile):
                cuemol_style(
                    "png",
                    filename=str(output / f"{name}.png"),
                    width=320,
                    height=240,
                    quiet=1,
                    _self=cmd,
                )
        a, b = (
            pixels((output / f"{name}.png").read_bytes()).astype(float)
            for name in ("single", "tiled")
        )
        error = float(np.abs(a - b).mean())
        assert error < 0.1, error
        from cuemol_style_in_pymol.sampling import camera

        cmd.set("orthoscopic", 0)
        cmd.turn("y", 23)
        pump()
        np.testing.assert_allclose(
            next(manager.active_drawings()).matrices[1], camera(cmd)[1], atol=2e-6
        )
        cmd.set("orthoscopic", 1)
        keys = (
            gl.GL_PIXEL_PACK_BUFFER_BINDING,
            gl.GL_PIXEL_UNPACK_BUFFER_BINDING,
            gl.GL_PACK_ROW_LENGTH,
            gl.GL_PACK_SKIP_ROWS,
            gl.GL_PACK_SKIP_PIXELS,
            gl.GL_UNPACK_ROW_LENGTH,
            gl.GL_UNPACK_SKIP_ROWS,
            gl.GL_UNPACK_SKIP_PIXELS,
            gl.GL_FRAMEBUFFER_BINDING,
            gl.GL_MATRIX_MODE,
            gl.GL_ACTIVE_TEXTURE,
        )
        with widget:
            gl.glPushClientAttrib(gl.GL_CLIENT_ALL_ATTRIB_BITS)
            old_pack = int(gl.glGetIntegerv(gl.GL_PIXEL_PACK_BUFFER_BINDING))
            old_unpack = int(gl.glGetIntegerv(gl.GL_PIXEL_UNPACK_BUFFER_BINDING))
            handles = gl.glGenBuffers(2)
            try:
                for target, handle in zip(
                    (gl.GL_PIXEL_PACK_BUFFER, gl.GL_PIXEL_UNPACK_BUFFER), handles
                ):
                    gl.glBindBuffer(target, int(handle))
                    gl.glBufferData(target, 64, None, gl.GL_STREAM_DRAW)
                for key in keys[2:8]:
                    gl.glPixelStorei(key, 7)
                before = [int(gl.glGetIntegerv(key)) for key in keys]
                matrix = np.asarray(gl.glGetDoublev(gl.GL_PROJECTION_MATRIX)).copy()
                manager.pool.draw(entry.drawings["protein"][0])
                assert before == [int(gl.glGetIntegerv(key)) for key in keys]
                np.testing.assert_array_equal(
                    matrix, gl.glGetDoublev(gl.GL_PROJECTION_MATRIX)
                )
                assert gl.glGetError() == gl.GL_NO_ERROR
            finally:
                gl.glBindBuffer(gl.GL_PIXEL_PACK_BUFFER, old_pack)
                gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, old_unpack)
                gl.glDeleteBuffers(2, handles)
                gl.glPopClientAttrib()
        cuemol_style("reset", quiet=1, _self=cmd)
        same_snapshot(cmd, baseline)
        for state in (2, 3):
            cmd.create("protein", "protein", 1, state)
        entry = cuemol_style("richardson", transparency=0.35, quiet=1, _self=cmd)
        pump()
        previous = next(manager.active_drawings()).sample_key
        cmd.turn("y", 47)
        manager.prepare_view()
        assert previous != next(manager.active_drawings()).sample_key
        assert cmd.get_setting_int("cgo_lighting", "cuemol_alpha_1") == 0
        for state in (3, 1, 2):
            cmd.frame(state)
            manager.prepare_view()
            assert entry.drawings["protein"][state - 1].sampled_bytes > 0
            assert all(
                not d.sampled_bytes
                for d in entry.drawings["protein"]
                if d.state != state
            )
            assert cmd.get_setting_int("cgo_lighting", "cuemol_alpha_1") == 0
        cuemol_style(
            "png",
            filename=str(output / "alpha.png"),
            width=320,
            height=240,
            quiet=1,
            _self=cmd,
        )
        saved = (output / "alpha.png").read_bytes()
        entry.options["cache_mb"] = 1
        try:
            cuemol_style(
                "png",
                filename=str(output / "alpha.png"),
                width=640,
                height=480,
                quiet=1,
                _self=cmd,
            )
        except CmdException:
            pass
        else:
            raise AssertionError("The sample budget must reject oversized exports")
        assert (output / "alpha.png").read_bytes() == saved
        assert not manager.pool.precise
        cuemol_style("reset", quiet=1, _self=cmd)
        cmd.delete("all")
        cmd.set("opaque_background", 1)
        cmd.set("ray_opaque_background", 1)
        for name, color, pos in (
            ("front", "red", [0, 0, 1]),
            ("rear", "blue", [0.5, 0, 0]),
        ):
            cmd.pseudoatom(name, pos=pos, elem="C", color=color)
            cmd.show_as("spheres", name)
        cmd.orient()
        cmd.zoom(buffer=3)
        cmd.clip("slab", 150)
        for name, transparency in (("front", 0.45), ("rear", 0.3)):
            cuemol_style(
                "cpk",
                name=name + "_style",
                selection=name,
                color="keep",
                transparency=transparency,
                quiet=1,
                _self=cmd,
            )
        pump()
        for mode in ("png", "ray"):
            cuemol_style(
                mode,
                filename=str(output / f"overlap-{mode}.png"),
                width=320,
                height=240,
                quiet=1,
                _self=cmd,
            )
        gpu, ray = (
            pixels((output / f"overlap-{mode}.png").read_bytes()).astype(float)
            for mode in ("png", "ray")
        )
        foreground = np.any(np.minimum(gpu, ray) < 250, axis=2)
        overlap_error = float(np.abs(gpu - ray)[foreground].mean())
        assert overlap_error < 1, overlap_error
        cuemol_style("reset", quiet=1, _self=cmd)
        cmd.delete("all")
        cmd.fab("AAAAAAAAAAAA", name="sheet", ss=2)
        cmd.remove("hydro")
        cmd.alter("sheet", 'ss="S"')
        cmd.show_as("cartoon")
        cmd.bg_color("black")
        cmd.orient()
        cmd.turn("x", 45)
        cmd.zoom(buffer=1)
        cmd.clip("slab", 150)
        cuemol_style("toon1", quiet=1, _self=cmd)
        for index in range(12):
            cmd.turn("y", 30)
            pump(0.03)
            cuemol_style(
                "png",
                filename=str(output / f"sheet-rotation-{index:02d}.png"),
                width=320,
                height=240,
                quiet=1,
                _self=cmd,
            )
            assert not next(manager.active_drawings()).error
        cuemol_style("reset", quiet=1, _self=cmd)
        pencil_error = check_live_pencil(cmd, widget, pump)
    report = {
        "tile_mean_error_255": error,
        "live_tile_mean_errors_255": live_errors,
        "live_pencil_native_mean_error_255": pencil_error,
        "large_perspective_first_draw": True,
        "transfer_state_restored": True,
        "alpha_camera_and_nonsequential_states": True,
        "failed_export_preserved_file": True,
        "overlap_gpu_ray_foreground_mae_255": overlap_error,
        "sheet_rotation_views": 12,
    }
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path(".cache/appearance-regression")
    )
    check(parser.parse_args().output)
