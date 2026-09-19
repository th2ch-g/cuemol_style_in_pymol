"""Native CGO interoperability and transactional image export."""

from contextlib import contextmanager
from io import BytesIO
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import uuid4

import numpy as np

from .materials import bake, srgb_decode, srgb_encode
from .mesh import unit


def downsample_png(data, width, height):
    """Average the 3x sample grid in premultiplied RGBA without interpolation."""
    from PIL import Image

    with Image.open(BytesIO(data)) as image:
        if image.size != (width * 3, height * 3):
            raise RuntimeError("PyMOL returned an unexpected supersampled image size")
        pixels = np.asarray(image.convert("RGBA"), dtype=np.float32)
        pixels[:, :, :3] *= pixels[:, :, 3:4] / 255
        pixels = pixels.reshape(height, 3, width, 3, 4).mean(axis=(1, 3))
        pixels[:, :, :3] *= np.divide(
            255,
            pixels[:, :, 3:4],
            out=np.zeros_like(pixels[:, :, 3:4]),
            where=pixels[:, :, 3:4] > 0,
        )
        result = Image.fromarray(np.uint8(np.clip(np.floor(pixels + 0.5), 0, 255)))
        output = BytesIO()
        options = {
            key: image.info[key] for key in ("dpi", "icc_profile") if key in image.info
        }
        result.save(output, format="PNG", **options)
        return output.getvalue()


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


def cgo_mesh(piece, material, rotation=None, sampled=None, *, normals=True):
    from pymol.cgo import ALPHA, BEGIN, COLOR, END, NORMAL, TRIANGLES, VERTEX

    mesh = piece.mesh if sampled is None else sampled[0]
    ids = mesh.faces.ravel()
    colors = bake(mesh, material, rotation) if sampled is None else sampled[1]
    values = np.empty((len(ids), 12 if normals else 8), np.float32)
    offset = 4 if normals else 0
    if normals:
        values[:, 0] = NORMAL
        values[:, 1:4] = mesh.normals[ids]
    values[:, offset] = COLOR
    values[:, offset + 1 : offset + 4] = colors[ids]
    values[:, offset + 4] = VERTEX
    values[:, offset + 5 : offset + 8] = mesh.vertices[ids]
    return [ALPHA, mesh.opacity, BEGIN, TRIANGLES, *values.ravel().tolist(), END]


def ray_proxy(drawing):
    """Retain opaque geometry for native ray without drawing it in OpenGL.

    PyMOL 3.1 consumes the public TRIANGLE opcode in its ray tracer and
    ignores it in both OpenGL CGO paths. Transparent bodies already have a
    native CGO, so including them here would render their opacity twice.
    Materials use molecular-space samples, independent of the active camera.
    The dedicated export adds camera-dependent contour geometry separately.
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
        size += sum(part.cgo_nbytes for part in drawing.native)
        opaque = [p for p in drawing.pieces if p.mesh.opacity >= 0.999999]
        if opaque:
            size += 4 * (2 + sum(28 * len(p.mesh.faces) for p in opaque))
        alpha = sum(
            4 * (5 + 36 * len(p.mesh.faces))
            for p in drawing.pieces
            if p.mesh.opacity < 0.999999
        )
        size += max(alpha, drawing.sampled_bytes)
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


def sampled_cgo(drawing, meshes, matrices, budget, opacity=1, native_fog=None):
    from .raster import colors, pixel_mesh
    from .sampling import cancel_native_fog

    image, depth = colors(drawing, meshes, matrices, budget)
    mesh = pixel_mesh(image, depth, matrices, budget, opacity, drawing.background)
    rgb = cancel_native_fog(mesh.colors, mesh, matrices, drawing.background, native_fog)
    return cgo_mesh(None, "nolighting", sampled=(mesh, rgb), normals=False)


def ray_cgo(drawing, cmd, width=0, height=0, budget=2048 * 1024**2):
    from .raster import layers
    from .sampling import camera

    matrices = camera(cmd, width, height, ray=True)
    drawing.background = tuple(cmd.get_color_tuple(cmd.get("bg_rgb")))
    view = cmd.get_view()
    drawing.fog = (-view[11], -view[11] + (view[16] - view[15]) / 2)
    result = []
    for opacity, meshes in layers(drawing).items():
        result.extend(sampled_cgo(drawing, meshes, matrices, budget, opacity))
        if len(result) * 4 > budget:
            raise ValueError("Camera-dependent ray CGO exceeds cache_mb")
    return result


def blend_images(base, passes):
    """Blend completed renderer passes in display space, including coverage."""
    from PIL import Image

    def pixels(data):
        with Image.open(BytesIO(data)) as image:
            values = np.asarray(image.convert("RGBA"), np.float32) / 255
        values[:, :, :3] = srgb_encode(values[:, :, :3])
        return values

    result = (1 - sum(alpha for alpha, _ in passes)) * pixels(base)
    for alpha, data in passes:
        result += alpha * pixels(data)
    result[:, :, :3] = srgb_decode(result[:, :, :3])
    output = BytesIO()
    Image.fromarray(np.uint8(np.clip(np.floor(result * 255 + 0.5), 0, 255))).save(
        output, format="PNG"
    )
    return output.getvalue()


def render_group_passes(manager, width, height, ray, temporary, disabled, prepared):
    from .gpu import Drawing
    from .presets import Profile
    from .raster import colors, layers
    from .sampling import camera

    cmd = manager.cmd
    matrices = camera(cmd, width, height, ray=ray)
    size = matrices[2][2:].astype(int)
    view = cmd.get_view()
    background = tuple(cmd.get_color_tuple(cmd.get("bg_rgb")))
    fog = (-view[11], -view[11] + (view[16] - view[15]) / 2)
    if ray:
        from .ray import can_compose, compose

        if can_compose(manager):
            drawings = list(manager.active_drawings())
            for drawing in drawings:
                drawing.background, drawing.fog = background, fog
            return compose(manager, drawings, matrices, disabled)
    groups = []
    callbacks = []
    for drawing in manager.active_drawings():
        drawing.background, drawing.fog = background, fog
        for opacity, meshes in layers(drawing).items():
            if not ray and opacity >= 0.999999:
                continue
            name = "_cuemol_export_" + uuid4().hex
            temporary.append(name)
            if ray:
                cmd.load_cgo(
                    sampled_cgo(drawing, meshes, matrices, drawing.sample_budget),
                    name,
                    zoom=0,
                )
                cmd.set("cgo_lighting", 0, name)
            else:
                image, depth = colors(drawing, meshes, matrices, drawing.sample_budget)
                projection = matrices[1]
                depth = (
                    0.5
                    + 0.5
                    * (-projection[2, 2] * depth + projection[2, 3])
                    / np.maximum(-projection[3, 2] * depth + projection[3, 3], 1e-8)
                ).astype(np.float32)
                depth[image[:, :, 3] == 0] = 1
                layer = Drawing(
                    [],
                    Profile(material="nolighting"),
                    drawing.edge_color,
                    manager.pool,
                    name=name,
                    background=background,
                    fog=(1e10, 2e10),
                    raster_image=image[::-1],
                    raster_depth=depth[::-1],
                    extent=drawing.extent,
                )
                callbacks.append(layer)
                prepared.append(layer)
                cmd.load_callback(layer, name, 1, 1, 0, 1, 0)
            if opacity < 0.999999:
                cmd.disable(name)
                groups.append((opacity, name))
    enabled = set(cmd.get_names("objects", enabled_only=1))
    for entry in manager.entries.values():
        alpha_names = {
            f"{entry.name}_alpha_{i}" for i in range(1, len(entry.drawings) + 1)
        }
        empty_shapes = {
            f"{entry.name}_shape_{i}"
            for i, drawings in enumerate(entry.drawings.values(), 1)
            if not any(p.mesh.opacity >= 0.999999 for d in drawings for p in d.pieces)
        }
        for name in entry.generated:
            if name in entry.native_objects:
                continue
            if name in enabled and (ray or name in alpha_names or name in empty_shapes):
                disabled.append(name)
                cmd.disable(name)
    manager.pool.raster_scale = 1

    def render():
        enabled = set(cmd.get_names("objects", enabled_only=1))
        active = [(d, d.raster_draws) for d in callbacks if d.name in enabled]
        (cmd.ray if ray else cmd.draw)(
            int(size[0]) * 3, int(size[1]) * 3, antialias=0, quiet=1
        )
        errors = [drawing.error for drawing in callbacks if drawing.error]
        if errors:
            raise RuntimeError(errors[0])
        if any(d.raster_draws == count for d, count in active):
            raise RuntimeError("The export layer does not match the active sample grid")
        return downsample_png(cmd.png(None, prior=1, quiet=1), *map(int, size))

    with settings(
        cmd,
        {
            "ray_trace_mode": 0,
            "ambient": 1,
            "direct": 0,
            "reflect": 0,
            "specular": 0,
            "light_count": 1,
            "ray_shadows": 0,
            "ray_trace_fog": 0,
            "ray_transparency_oblique": 0,
            "ray_transparency_contrast": 1,
            "ray_legacy_lighting": 0,
        }
        if ray
        else {},
    ):
        base = render()
        passes = []
        for opacity, name in groups:
            cmd.enable(name)
            try:
                passes.append((opacity, render()))
            finally:
                cmd.disable(name)
        return blend_images(base, passes) if passes else base


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
    raster_scale = manager.pool.raster_scale
    precise = manager.pool.precise
    manager.pool.precise = True
    manager.pool.show_selection = False
    temporary, disabled, prepared = [], [], []
    busy = manager.busy
    try:
        if playing:
            cmd.mstop()
        if not ray and manager.widget is None:
            raise ValueError(
                "GPU PNG export requires the PyMOL Qt GUI; use cuemol_style ray in headless mode"
            )
        manager.prepare_view(width, height, native=False)
        manager.busy = True
        alpha = any(
            p.mesh.opacity < 0.999999
            for d in manager.active_drawings()
            for p in d.pieces
        )
        if ray or alpha:
            data = render_group_passes(
                manager, width, height, ray, temporary, disabled, prepared
            )
        else:
            cmd.draw(width, height, antialias=0, quiet=1)
            data = cmd.png(None, prior=1, quiet=1)
        if not ray:
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
        manager.busy = busy
        manager.pool.show_selection = show_selection
        manager.pool.raster_scale = raster_scale
        manager.pool.precise = precise
        for name in temporary:
            cmd.delete(name)
        if prepared:
            manager.pool.hatch.key = None
        for name in disabled:
            cmd.enable(name)
        if cmd.get_frame() != frame:
            cmd.frame(frame)
        if sculpting:
            cmd.set("sculpting", sculpting)
        if playing:
            cmd.mplay()
        manager.prepare_view()
        cmd.refresh()
    return str(path)
