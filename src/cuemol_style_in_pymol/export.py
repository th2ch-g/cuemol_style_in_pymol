"""Native CGO interoperability and transactional image export."""

from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import uuid4

import numpy as np

from .materials import bake
from .mesh import unit


@contextmanager
def settings(cmd, values):
    saved = {key: cmd.get_setting_tuple(key)[1] for key in values}
    try:
        for key, value in values.items():
            cmd.set(key, value)
        yield
    finally:
        for key, value in saved.items():
            cmd.set(key, value if len(value) > 1 else value[0])


def cgo_mesh(piece, material, rotation=None):
    from pymol.cgo import ALPHA, BEGIN, COLOR, END, NORMAL, TRIANGLES, VERTEX

    mesh = piece.mesh
    ids = mesh.faces.ravel()
    colors = bake(mesh, material, rotation)
    values = np.empty((len(ids), 12), float)
    values[:, 0] = NORMAL
    values[:, 1:4] = mesh.normals[ids]
    values[:, 4] = COLOR
    values[:, 5:8] = colors[ids]
    values[:, 8] = VERTEX
    values[:, 9:12] = mesh.vertices[ids]
    return [ALPHA, mesh.opacity, BEGIN, TRIANGLES, *values.ravel().tolist(), END]


def ray_proxy(drawing):
    """Retain opaque geometry for native ray without drawing it in OpenGL.

    PyMOL 3.1 consumes the public TRIANGLE opcode in its ray tracer and
    ignores it in both OpenGL CGO paths. Transparent bodies already have a
    native CGO, so including them here would render their opacity twice.
    Materials use molecular-space samples, independent of the active camera.
    The dedicated export adds camera-dependent outline cylinders separately.
    """
    from pymol.cgo import ALPHA, TRIANGLE

    parts = [np.asarray([ALPHA, 1.0], np.float32)]
    for piece in drawing.pieces:
        mesh = piece.mesh
        if mesh.opacity < 0.999999:
            continue
        # Match the vertex order used by PyMOL's BEGIN/TRIANGLES ray path.
        faces = mesh.faces[:, ::-1]
        values = np.empty((len(mesh.faces), 28), np.float32)
        values[:, 0] = TRIANGLE
        values[:, 1:10] = mesh.vertices[faces].reshape(-1, 9)
        values[:, 10:19] = mesh.normals[faces].reshape(-1, 9)
        colors = bake(mesh, drawing.profile.material)
        values[:, 19:28] = colors[faces].reshape(-1, 9)
        parts.append(values.ravel())
    return np.concatenate(parts).tolist()


def native_cgo_bytes(drawings):
    """Estimate float32 command payloads, excluding native allocator overhead."""
    size = 0
    for drawing in drawings:
        opaque = [p for p in drawing.pieces if p.mesh.opacity >= 0.999999]
        if opaque:
            size += 4 * (2 + sum(28 * len(p.mesh.faces) for p in opaque))
        for piece in drawing.pieces:
            if piece.mesh.opacity < 0.999999:
                size += 4 * (5 + 36 * len(piece.mesh.faces))
    return size


def view_matrix(cmd):
    view = np.asarray(cmd.get_view())
    matrix = np.eye(4)
    matrix[:3, :3] = view[:9].reshape(3, 3).T
    matrix[:3, 3] = view[9:12] - matrix[:3, :3] @ view[12:15]
    return matrix


def visible_edges(edges, modelview, orthoscopic, creases):
    if not len(edges):
        return edges
    rotation = modelview[:3, :3]
    mid = 0.5 * (edges[:, :3] + edges[:, 3:6]) @ rotation.T + modelview[:3, 3]
    direction = np.tile([0.0, 0.0, 1.0], (len(mid), 1)) if orthoscopic else unit(-mid)
    a, b = edges[:, 6:9] @ rotation.T, edges[:, 9:12] @ rotation.T
    fa, fb = np.sum(a * direction, axis=1), np.sum(b * direction, axis=1)
    keep = fa * fb <= 0
    if creases:
        keep |= (np.sum(a * b, axis=1) < 0.5) & (np.maximum(fa, fb) > 0)
    return edges[keep]


def ray_cgo(drawing, cmd):
    from pymol.cgo import ALPHA, CYLINDER

    matrix = view_matrix(cmd)
    result = []
    for piece in drawing.pieces:
        result.extend(cgo_mesh(piece, drawing.profile.material, matrix[:3, :3]))
        if drawing.profile.edges != "none":
            edges = visible_edges(
                piece.edges,
                matrix,
                cmd.get_setting_int("orthoscopic"),
                drawing.profile.edges == "edges",
            )
            values = np.empty((len(edges), 14), float)
            values[:, 0] = CYLINDER
            values[:, 1:7] = edges[:, :6]
            values[:, 7] = drawing.profile.edge_width / 2
            values[:, 8:11] = drawing.edge_color
            values[:, 11:14] = drawing.edge_color
            result.extend([ALPHA, 1.0, *values.ravel().tolist()])
    return result


def image(manager, filename, width, height, ray):
    if not filename:
        raise ValueError("filename is required for PNG and ray export")
    width, height = int(width), int(height)
    if width < 0 or height < 0 or width > 32768 or height > 32768:
        raise ValueError("Image dimensions must be between 0 and 32768 pixels")
    path = Path(filename).expanduser()
    if path.suffix.lower() != ".png":
        path = Path(str(path) + ".png")
    if not path.parent.is_dir():
        raise ValueError("The output directory does not exist")
    cmd = manager.cmd
    playing = cmd.get_movie_playing()
    sculpting = cmd.get_setting_int("sculpting")
    frame = cmd.get_frame()
    show_selection = manager.pool.show_selection
    manager.pool.show_selection = False
    temporary, disabled = [], []
    try:
        if playing:
            cmd.mstop()
        if ray:
            for drawing in manager.active_drawings():
                name = "_cuemol_ray_" + uuid4().hex
                temporary.append(name)
                cmd.load_cgo(ray_cgo(drawing, cmd), name, zoom=0)
            for entry in manager.entries.values():
                for name in entry.generated:
                    if name in cmd.get_names("objects", enabled_only=1):
                        disabled.append(name)
                        cmd.disable(name)
            # Edges are explicit geometry, avoiding an extra global outline pass.
            with settings(cmd, {"ray_trace_mode": 0}):
                data = cmd.png(None, width, height, ray=1, quiet=1)
        else:
            if manager.widget is None:
                raise ValueError(
                    "GPU PNG export requires the PyMOL Qt GUI; use cuemol_style ray in headless mode"
                )
            cmd.draw(width, height, quiet=1)
            data = cmd.png(None, prior=1, quiet=1)
            errors = [
                d.error
                for e in manager.entries.values()
                for ds in e.drawings.values()
                for d in ds
                if d.error
            ]
            if errors:
                raise RuntimeError(errors[0])
        if not isinstance(data, bytes) or not data.startswith(b"\x89PNG\r\n\x1a\n"):
            raise RuntimeError("PyMOL did not return a PNG image")
        with NamedTemporaryFile(dir=path.parent, suffix=".png", delete=False) as stream:
            scratch = Path(stream.name)
            try:
                stream.write(data)
                stream.flush()
                scratch.replace(path)
            finally:
                scratch.unlink(missing_ok=True)
    finally:
        manager.pool.show_selection = show_selection
        for name in temporary:
            cmd.delete(name)
        for name in disabled:
            cmd.enable(name)
        if cmd.get_frame() != frame:
            cmd.frame(frame)
        if sculpting:
            cmd.set("sculpting", sculpting)
        if playing:
            cmd.mplay()
        cmd.rebuild()
        cmd.refresh()
    return str(path)
