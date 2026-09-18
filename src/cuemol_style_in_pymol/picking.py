"""Qt picking adapter which leaves native mouse gestures intact."""

import numpy as np

from .mesh import ray_hits, unit
from .source import atom_selection


def hit_at(drawings, x, y, depth=1.0):
    """Intersect screen rays with the same indexed geometry used for drawing."""
    hits = []
    for drawing in drawings:
        if drawing.matrices is None:
            continue
        modelview, projection, viewport = drawing.matrices
        vx, vy, width, height = viewport
        if not (vx <= x < vx + width and vy <= y < vy + height):
            continue
        matrix = projection @ modelview
        inverse = np.linalg.inv(matrix)
        ndc = [2 * (x - vx) / width - 1, 2 * (y - vy) / height - 1]
        near = inverse @ [*ndc, -1, 1]
        far = inverse @ [*ndc, 1, 1]
        origin = near[:3] / near[3]
        direction = unit(far[:3] / far[3] - origin)
        for piece in drawing.pieces:
            hit = ray_hits(piece.mesh, origin, direction)
            if hit is None:
                continue
            distance, owner = hit
            point = origin + direction * distance
            clip = matrix @ [*point, 1]
            z = (clip[2] / clip[3] + 1) / 2
            # Native opaque geometry must occlude custom picking as well.
            if 0 <= z <= min(1.0, depth + 2e-5):
                hits.append((z, piece.atoms[owner]))
    return min(hits, key=lambda value: value[0])[1] if hits else None


def select_atom(cmd, atom, additive=False):
    modes = {
        0: "",
        1: "byres",
        2: "bychain",
        3: "byseg",
        4: "byobject",
        5: "bymolecule",
        6: "bycalpha",
    }
    mode = modes.get(cmd.get_setting_int("mouse_selection_mode"), "")
    with atom_selection(cmd, [(atom.model, atom.index)]) as sele:
        expression = f"{mode} ({sele})"
        if additive and "sele" in cmd.get_names("selections"):
            expression = f"(sele) or ({expression})"
        cmd.select("sele", expression, quiet=1)
    cmd.enable("sele")
    cmd.refresh()


def attach(manager):
    try:
        from pymol.Qt import QtCore, QtGui, QtWidgets
    except ImportError:
        return None, None
    app = QtWidgets.QApplication.instance()
    if app is None:
        return None, None
    widgets = [
        w
        for w in app.allWidgets()
        if hasattr(w, "makeCurrent")
        and getattr(getattr(w, "cmd", None), "_COb", None) == manager.cmd._COb
    ]
    if not widgets:
        return None, None
    widget = widgets[0]

    class Events(QtCore.QObject):
        def __init__(self):
            super().__init__(widget)
            self.press = None
            self.forwarding = False
            self.closed = False
            self.warned = False
            self.timer = QtCore.QTimer(self)
            self.timer.timeout.connect(self.maintain)
            self.timer.start(100)
            self.settle = QtCore.QTimer(self)
            self.settle.setSingleShot(True)
            self.settle.timeout.connect(self.refine)
            handlers = getattr(widget, "_pymol_style_pickers", [])
            widget._pymol_style_pickers = [*handlers, self]
            widget.installEventFilter(self)

        def interact(self):
            manager.pool.interacting = True
            self.settle.start(150)

        def refine(self):
            if self.closed or not manager.pool.interacting:
                return
            if manager.busy:
                self.settle.start(150)
                return
            manager.pool.interacting = False
            # Invalidate PyMOL's cached scene outside its drawing callback.
            manager.cmd.refresh()
            widget.update()

        def maintain(self):
            try:
                manager.maintenance()
            except Exception as exc:  # noqa: BLE001 - GUI callback boundary.
                self.timer.stop()
                print(
                    f" cuemol_style: view maintenance stopped: {exc}; use refresh or reset."
                )

        def close(self):
            self.closed = True
            self.timer.stop()
            self.settle.stop()
            manager.pool.interacting = False
            widget._pymol_style_pickers = [
                handler
                for handler in widget._pymol_style_pickers
                if handler is not self
            ]
            widget.removeEventFilter(self)
            self.deleteLater()

        def depth(self, x, y):
            from OpenGL import GL as gl
            from OpenGL.GL.EXT.framebuffer_object import glBindFramebufferEXT

            with widget:
                old = int(gl.glGetIntegerv(gl.GL_FRAMEBUFFER_BINDING))
                bind = (
                    gl.glBindFramebuffer
                    if gl.glBindFramebuffer
                    else glBindFramebufferEXT
                )
                try:
                    bind(
                        gl.GL_FRAMEBUFFER,
                        widget.defaultFramebufferObject()
                        if hasattr(widget, "defaultFramebufferObject")
                        else 0,
                    )
                    value = gl.glReadPixels(
                        int(x), int(y), 1, 1, gl.GL_DEPTH_COMPONENT, gl.GL_FLOAT
                    )
                    return float(np.asarray(value).ravel()[0])
                finally:
                    bind(gl.GL_FRAMEBUFFER, old)

        def hit(self, event):
            ratio = widget.devicePixelRatioF()
            x = int(event.pos().x() * ratio) + 0.5
            y = int((widget.height() - event.pos().y()) * ratio) + 0.5
            drawings = list(manager.active_drawings())
            depth = self.depth(x, y)
            atom = hit_at(drawings, x, y, depth)
            if atom is not None:
                return atom
            # Contour compositing retains the nearest of its nine subpixels.
            # Test that footprint before treating the depth as a native occluder.
            for dy in (-1 / 3, 0, 1 / 3):
                for dx in (-1 / 3, 0, 1 / 3):
                    if dx or dy:
                        atom = hit_at(drawings, x + dx, y + dy, depth)
                        if atom is not None:
                            return atom
            return None

        def pick(self, atom, additive):
            if self.closed or manager.busy:
                return
            try:
                select_atom(manager.cmd, atom, additive)
                manager.update_selection()
            except Exception as exc:  # noqa: BLE001 - GUI callback boundary.
                if not self.warned:
                    print(
                        f" cuemol_style: picking failed: {exc}; command-line selections remain available."
                    )
                    self.warned = True

        def forward(self, press):
            # Replayed native presses must bypass every sibling event filter.
            widget._pymol_style_forwarding = True
            self.forwarding = True
            try:
                QtWidgets.QApplication.sendEvent(widget, press)
            finally:
                self.forwarding = False
                widget._pymol_style_forwarding = False

        def eventFilter(self, watched, event):
            if (
                self.forwarding
                or getattr(widget, "_pymol_style_forwarding", False)
                or self.closed
                or manager.busy
            ):
                return False
            kind = event.type()
            if kind in (QtCore.QEvent.FocusOut, QtCore.QEvent.Hide):
                self.press = None
            if kind == QtCore.QEvent.Wheel or (
                kind == QtCore.QEvent.MouseMove and event.buttons()
            ):
                self.interact()
            elif (
                kind
                in (
                    QtCore.QEvent.MouseButtonRelease,
                    QtCore.QEvent.FocusOut,
                    QtCore.QEvent.Hide,
                )
                and manager.pool.interacting
            ):
                self.settle.start(150)
            if event.type() == QtCore.QEvent.Paint:
                try:
                    manager.prepare_view()
                except Exception as exc:  # noqa: BLE001 - GUI callback boundary.
                    for drawing in manager.active_drawings():
                        drawing.error = str(exc)
                return False
            if event.type() in (
                QtCore.QEvent.MouseButtonPress,
                QtCore.QEvent.MouseButtonDblClick,
            ):
                if event.button() == QtCore.Qt.LeftButton and not (
                    event.modifiers()
                    & (QtCore.Qt.ControlModifier | QtCore.Qt.AltModifier)
                ):
                    # Wait until release to distinguish clicks from drags. A
                    # drag must not wait for GPU depth readback or ray tests.
                    self.press = QtGui.QMouseEvent(event)
                    return True
            elif event.type() == QtCore.QEvent.MouseMove and self.press is not None:
                press = self.press
                if (event.pos() - press.pos()).manhattanLength() <= 4:
                    return True
                self.press = None
                self.forward(press)
            elif (
                event.type() == QtCore.QEvent.MouseButtonRelease
                and self.press is not None
                and event.button() == QtCore.Qt.LeftButton
            ):
                press = self.press
                self.press = None
                # A sibling style may own the visible custom object. Native
                # depth rejects covered geometry before choosing its picker.
                for handler in [
                    self,
                    *(h for h in widget._pymol_style_pickers if h is not self),
                ]:
                    try:
                        atom = handler.hit(press)
                    except Exception:  # noqa: BLE001 - GUI callback boundary.
                        continue
                    if atom is not None:
                        handler.pick(
                            atom, bool(event.modifiers() & QtCore.Qt.ShiftModifier)
                        )
                        return True
                # Let native geometry and background clicks keep PyMOL's normal
                # selection behavior; custom clicks never reach its picker.
                self.forward(press)
            return False

    return widget, Events()
