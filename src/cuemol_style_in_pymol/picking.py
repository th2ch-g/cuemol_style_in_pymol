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
            widget.installEventFilter(self)

        def maintain(self):
            try:
                manager.maintenance()
            except Exception as exc:
                self.timer.stop()
                print(
                    f" cuemol_style: view maintenance stopped: {exc}; use refresh or reset."
                )

        def close(self):
            self.closed = True
            self.timer.stop()
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
            x = event.pos().x() * ratio
            y = (widget.height() - event.pos().y()) * ratio
            return hit_at(list(manager.active_drawings()), x, y, self.depth(x, y))

        def pick(self, atom, additive):
            if self.closed or manager.busy:
                return
            try:
                select_atom(manager.cmd, atom, additive)
                manager.update_selection()
            except Exception as exc:
                if not self.warned:
                    print(
                        f" cuemol_style: picking failed: {exc}; command-line selections remain available."
                    )
                    self.warned = True

        def eventFilter(self, watched, event):
            if self.forwarding or self.closed or manager.busy:
                return False
            if event.type() in (
                QtCore.QEvent.MouseButtonPress,
                QtCore.QEvent.MouseButtonDblClick,
            ):
                if event.button() == QtCore.Qt.LeftButton and not (
                    event.modifiers()
                    & (QtCore.Qt.ControlModifier | QtCore.Qt.AltModifier)
                ):
                    try:
                        atom = self.hit(event)
                    except Exception:
                        return False
                    if atom is not None:
                        # Hold only custom presses. Replaying a native click
                        # would trigger PyMOL's deferred background deselect.
                        self.press = QtGui.QMouseEvent(event), atom
                        return True
            elif event.type() == QtCore.QEvent.MouseMove and self.press is not None:
                press, _ = self.press
                if (event.pos() - press.pos()).manhattanLength() <= 4:
                    return True
                self.press = None
                self.forwarding = True
                try:
                    QtWidgets.QApplication.sendEvent(widget, press)
                finally:
                    self.forwarding = False
            elif event.type() == QtCore.QEvent.MouseButtonRelease:
                if self.press is not None and event.button() == QtCore.Qt.LeftButton:
                    _, atom = self.press
                    self.press = None
                    self.pick(atom, bool(event.modifiers() & QtCore.Qt.ShiftModifier))
                    return True
            return False

    return widget, Events()
