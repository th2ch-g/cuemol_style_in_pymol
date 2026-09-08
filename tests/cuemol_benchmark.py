"""Deterministic real-GUI performance checks for cuemol_style."""

from copy import deepcopy
import json
import resource
import sys
from time import perf_counter
from unittest.mock import patch

import numpy as np


def benchmark(cmd, widget, pump, output, style, report, structure=None):
    from chempy.models import Indexed
    from OpenGL import GL as gl
    from pymol.Qt import QtWidgets
    from cuemol_style_in_pymol.controller import manager_for

    print("Preparing the 500-residue / 100-state benchmark", flush=True)
    if structure:
        cmd.load(str(structure), "peptide")
    else:
        cmd.fab("ACDEFGHIKLMN", name="peptide", ss=1)
        cmd.alter("peptide and resi 8-10", 'ss="S"')
    template = cmd.get_model("peptide and polymer.protein and name N+CA+C+O")
    residues = list(dict.fromkeys((a.chain, a.resi) for a in template.atom))
    if not residues:
        raise ValueError("The benchmark requires a protein backbone")
    lookup = {key: i for i, key in enumerate(residues)}
    base_center = np.asarray(template.get_coord_list()).mean(axis=0)
    model = Indexed()
    total, block = 0, 0
    while total < 500:
        take = min(len(residues), 500 - total)
        included = {
            i: len(model.atom) + j
            for j, i in enumerate(
                i
                for i, a in enumerate(template.atom)
                if lookup[(a.chain, a.resi)] < take
            )
        }
        for i in included:
            atom = deepcopy(template.atom[i])
            rank = lookup[(atom.chain, atom.resi)]
            atom.resi = str(total + rank + 1)
            atom.chain = chr(65 + block % 26)
            atom.segi = str(block // 26)
            atom.coord = (
                np.asarray(atom.coord)
                - base_center
                + [block % 4 * 35, block // 4 * 35, 0]
            ).tolist()
            model.atom.append(atom)
        for original in template.bond:
            if all(i in included for i in original.index):
                bond = deepcopy(original)
                bond.index = [included[i] for i in original.index]
                model.bond.append(bond)
        total += take
        block += 1
    cmd.delete("all")
    cmd.load_model(model, "trajectory")
    base = cmd.get_coords("trajectory")
    phase = np.arange(len(base))[:, None] * 0.02
    for state in range(2, 101):
        coords = base + 0.1 * np.sin(phase + state * 0.12) * np.array([[1.0, 0.5, 0.3]])
        cmd.load_coordset(coords, "trajectory", state)
    assert cmd.count_atoms("trajectory and name CA") == 500
    assert cmd.count_states("trajectory") == 100
    cmd.show_as("cartoon")
    cmd.deselect()
    cmd.reset()
    cmd.zoom("trajectory", 3)
    ratio = widget.devicePixelRatioF()
    cmd.set("internal_gui", 0)
    cmd.set("internal_feedback", 0)
    cmd.set("seq_view", 0)
    cmd.set("movie_panel", 0)
    cmd.set("internal_prompt", 0)
    widget.setFixedSize(round(1280 / ratio), round(720 / ratio))
    widget.window().adjustSize()
    pump(0.2)
    # Recompute PyMOL's scene rectangle after hiding its internal panels.
    with widget:
        widget.resizeGL(widget.width(), widget.height())
    entry = style("richardson", quality="medium", quiet=1, _self=cmd)
    manager = manager_for(cmd)
    for _ in range(5):
        # Invalidate PyMOL's cached image after Qt has resized the widget.
        cmd.turn("y", 0)
        widget.repaint()
        pump(0.2)
        viewport = next(manager.active_drawings()).matrices[2]
        if tuple(viewport[2:]) == (1280, 720):
            break
    assert tuple(viewport[2:]) == (1280, 720), viewport
    print(
        f"Prepared {entry.nbytes / 1024**2:.1f} MiB of geometry in {entry.seconds:.2f} seconds",
        flush=True,
    )

    def render():
        widget.repaint()
        QtWidgets.QApplication.processEvents()
        with widget:
            gl.glFinish()

    # Warm all states before the timed passes, retaining the bounded GPU cache.
    start = perf_counter()
    for state in range(1, 101):
        cmd.frame(state)
        render()
        assert entry.drawings["trajectory"][state - 1].draws > 0
    warmup = perf_counter() - start
    cmd.frame(1)
    start = perf_counter()
    for _ in range(120):
        cmd.turn("y", 3)
        render()
    rotation_fps = 120 / (perf_counter() - start)
    start = perf_counter()
    for state in list(range(1, 101)) * 2:
        cmd.frame(state)
        render()
    playback_fps = 200 / (perf_counter() - start)
    cmd.mset("1 -100")
    cmd.set("movie_fps", 30)
    cmd.frame(1)
    seen = []
    original = manager.pool.draw

    def record(drawing):
        original(drawing)
        if not seen or seen[-1] != drawing.state:
            seen.append(drawing.state)

    with patch.object(manager.pool, "draw", record):
        cmd.mplay()
        start = perf_counter()
        pump(3)
        movie_seconds = perf_counter() - start
        cmd.mstop()
    stats = {
        "residues": 500,
        "states": 100,
        "quality": "medium",
        "viewport": [1280, 720],
        "input": "repeated peptide coordinates with deterministic synthetic state displacements",
        "prepare_seconds": entry.seconds,
        "warmup_seconds": warmup,
        "mesh_cache_mib": entry.nbytes / 1024**2,
        "native_cgo_mib": entry.cgo_nbytes / 1024**2,
        "gpu_cache_mib": manager.pool.bytes / 1024**2,
        "rotation_fps": rotation_fps,
        "state_switch_fps": playback_fps,
        "movie_drawn_states_per_second": len(seen) / movie_seconds,
    }
    print(json.dumps(stats, indent=2), flush=True)
    cmd.frame(1)
    style(
        "png", filename=str(output / "benchmark.png"), width=1280, height=720, _self=cmd
    )
    style("reset", _self=cmd)
    cmd.mset("")
    cmd.create("surface_case", "trajectory", 1, 1)
    cmd.delete("trajectory")
    cmd.frame(1)
    surface = style(
        "matte", representation="surface", quality="medium", quiet=1, _self=cmd
    )
    pump(0.1)
    stats["surface_prepare_seconds"] = surface.seconds
    stats["surface_mesh_mib"] = surface.nbytes / 1024**2
    stats["peak_process_rss_mib"] = resource.getrusage(
        resource.RUSAGE_SELF
    ).ru_maxrss / (1024**2 if sys.platform == "darwin" else 1024)
    style("reset", _self=cmd)
    cmd.delete("all")
    report["benchmark"] = stats
