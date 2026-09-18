"""Measure completed native mouse rotations in an isolated PyMOL window."""

import argparse
import importlib
import json
from pathlib import Path
from time import perf_counter

from check_cuemol_style import qt_session
from OpenGL import GL as gl
from pymol.Qt import QtCore, QtGui, QtWidgets

parser = argparse.ArgumentParser()
parser.add_argument("package", choices=["cuemol", "molstar", "chimerax"])
parser.add_argument("styles", nargs="+")
parser.add_argument("--output", type=Path, default=Path(".cache/drag.json"))
parser.add_argument("--structure", type=Path, required=True)
parser.add_argument("--width", type=int, default=3200)
parser.add_argument("--height", type=int, default=1800)
args = parser.parse_args()
module = importlib.import_module(args.package + "_style_in_pymol")
style = getattr(module, args.package + "_style")
report = {}
with qt_session() as (cmd, widget, pump):
    for key in (
        "internal_gui",
        "internal_feedback",
        "internal_prompt",
        "seq_view",
        "movie_panel",
    ):
        cmd.set(key, 0)
    ratio = widget.devicePixelRatioF()
    widget.setFixedSize(round(args.width / ratio), round(args.height / ratio))
    pump()
    with widget:
        widget.resizeGL(widget.width(), widget.height())
    cmd.load(str(args.structure), "molecule")
    cmd.show_as("cartoon")
    cmd.orient()
    cmd.zoom(buffer=3)
    cmd.set("orthoscopic", 0)
    view = cmd.get_view()
    for name in args.styles:
        cmd.set_view(view)
        t = perf_counter()
        entry = style(name, quiet=1, _self=cmd)
        prepare = perf_counter() - t
        pump(0.2)
        manager = vars(cmd._pymol).get("_" + args.package + "_style_manager")
        start = QtCore.QPointF(widget.width() / 2, widget.height() / 2)

        def event(kind, pos, button, buttons):
            QtWidgets.QApplication.sendEvent(
                widget,
                QtGui.QMouseEvent(kind, pos, button, buttons, QtCore.Qt.NoModifier),
            )

        t = perf_counter()
        event(
            QtCore.QEvent.MouseButtonPress,
            start,
            QtCore.Qt.LeftButton,
            QtCore.Qt.LeftButton,
        )
        press_ms = (perf_counter() - t) * 1000
        for warmup in range(3):
            pos = start + QtCore.QPointF(10 + warmup * 10, 10 + warmup * 6)
            event(
                QtCore.QEvent.MouseMove, pos, QtCore.Qt.NoButton, QtCore.Qt.LeftButton
            )
            widget.pymol.idle()
            with widget:
                widget.paintGL()
                gl.glFinish()
        total = perf_counter()
        times = []
        changed_frames = 0
        last_view = cmd.get_view()
        for i in range(40):
            t = perf_counter()
            pos = start + QtCore.QPointF(40 + 3 * i, 30 + 2 * i)
            event(
                QtCore.QEvent.MouseMove, pos, QtCore.Qt.NoButton, QtCore.Qt.LeftButton
            )
            widget.pymol.idle()
            this_view = cmd.get_view()
            changed_frames += this_view != last_view
            last_view = this_view
            with widget:
                widget.paintGL()
                gl.glFinish()
            times.append(perf_counter() - t)
        elapsed = perf_counter() - total
        event(
            QtCore.QEvent.MouseButtonRelease,
            pos,
            QtCore.Qt.LeftButton,
            QtCore.Qt.NoButton,
        )
        pump(0.4)
        data = {
            "prepare_seconds": prepare,
            "press_ms": press_ms,
            "completed_drag_frames_per_second": changed_frames / elapsed,
            "max_drag_frame_ms": max(times) * 1000,
            "view_changed": view != cmd.get_view(),
        }
        if manager:
            data["draw_errors"] = [
                d.error for d in manager.active_drawings() if d.error
            ]
            data["settled"] = not getattr(manager.pool, "interacting", False)
            data["mesh_mib"] = getattr(entry, "nbytes", 0) / 1024**2
            data["gpu_mib"] = manager.pool.bytes / 1024**2
            if args.package == "cuemol" and hasattr(manager.pool, "preview"):
                data["preview_mib"] = manager.pool.preview.bytes / 1024**2
                data["static_mib"] = manager.pool.live.bytes / 1024**2
            assert not data["draw_errors"] and data["settled"]
        data["changed_frames"] = changed_frames
        assert changed_frames >= 39, data
        report[name] = data
        print(name, data, flush=True)
        style("reset", quiet=1, _self=cmd)
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report, indent=2))
