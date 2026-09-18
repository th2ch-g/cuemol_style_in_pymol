"""Check native mouse gestures when multiple independent style filters coexist."""

from contextlib import ExitStack
from unittest.mock import patch

import numpy as np
from check_cuemol_style import qt_session
from pymol.Qt import QtCore, QtGui, QtWidgets
from PyQt5.QtTest import QTest

from cuemol_style_in_pymol import cuemol_style
from cuemol_style_in_pymol.controller import manager_for as cue_manager
from molstar_style_in_pymol import molstar_style
from molstar_style_in_pymol.controller import manager_for as mol_manager


def main():
    with qt_session() as (cmd, widget, pump):
        for reverse in (False, True):
            cmd.delete("all")
            cmd.pseudoatom("cue_atom", pos=[-3, 0, 0], elem="O", vdw=1.2)
            cmd.pseudoatom("mol_atom", pos=[3, 0, 0], elem="N", vdw=1.2)
            cmd.pseudoatom("native_atom", pos=[0, 0, 3], elem="C", vdw=1.2)
            cmd.show_as("spheres")
            cmd.reset()
            cmd.orient()
            cmd.zoom(buffer=3)
            cmd.set("mouse_selection_mode", 0)
            commands = [
                (cuemol_style, "cpk", "cue_atom"),
                (molstar_style, "spacefill", "mol_atom"),
            ]
            for command, profile, selection in (
                reversed(commands) if reverse else commands
            ):
                command(profile, selection=selection, quiet=1, _self=cmd)
            managers = [cue_manager(cmd), mol_manager(cmd)]
            pump(0.3)

            def point(world):
                drawing = next(managers[0].active_drawings())
                mv, projection, viewport = drawing.matrices
                clip = projection @ mv @ [*world, 1]
                ndc = clip[:3] / clip[3]
                ratio = widget.devicePixelRatioF()
                x = viewport[0] + (ndc[0] + 1) * viewport[2] / 2
                y = viewport[1] + (ndc[1] + 1) * viewport[3] / 2
                return QtCore.QPoint(
                    round(x / ratio), round(widget.height() - y / ratio)
                )

            for name in ("cue_atom", "mol_atom", "native_atom"):
                cmd.deselect()
                QTest.mouseClick(
                    widget, QtCore.Qt.LeftButton, pos=point(cmd.get_coords(name)[0])
                )
                pump(0.3)
                assert cmd.index("sele") == [(name, 1)], (
                    reverse,
                    name,
                    cmd.index("sele"),
                )
            cmd.deselect()
            before = cmd.get_view()
            start = point(cmd.get_coords("cue_atom")[0])
            finish = start + QtCore.QPoint(80, 40)
            with ExitStack() as stack:
                spies = [
                    stack.enter_context(patch.object(m.events, "hit")) for m in managers
                ]
                QTest.mousePress(widget, QtCore.Qt.LeftButton, pos=start)
                QtWidgets.QApplication.sendEvent(
                    widget,
                    QtGui.QMouseEvent(
                        QtCore.QEvent.MouseMove,
                        QtCore.QPointF(finish),
                        QtCore.Qt.NoButton,
                        QtCore.Qt.LeftButton,
                        QtCore.Qt.NoModifier,
                    ),
                )
                QTest.mouseRelease(widget, QtCore.Qt.LeftButton, pos=finish)
                for spy in spies:
                    spy.assert_not_called()
            pump(0.3)
            assert not np.array_equal(before, cmd.get_view())
            assert all(not m.pool.interacting for m in managers)
            for command, _, _ in commands:
                command("reset", name="all", quiet=1, _self=cmd)
            assert not widget._pymol_style_pickers
            print(
                f"Sibling picking, native picking, drag, and reset passed (reverse={reverse})"
            )


if __name__ == "__main__":
    main()
