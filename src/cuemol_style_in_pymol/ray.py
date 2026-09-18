"""Compose sampled ray layers without expanding pixels into native triangles."""

from io import BytesIO

import numpy as np
from PIL import Image


def can_compose(manager):
    """Keep native intersection and color conversion for mixed or gamma scenes."""
    cmd = manager.cmd
    if abs(cmd.get_setting_float("gamma") - 1) > 1e-6:
        return False
    generated = {name for e in manager.entries.values() for name in e.generated}
    for name in cmd.get_names("objects", enabled_only=1):
        if name in generated:
            continue
        kind = cmd.get_type(name)
        if kind == "object:molecule":
            visible = []
            cmd.iterate("%" + name, "out.append(reps)", space={"out": visible})
            if any(visible):
                return False
        elif kind not in ("object:group", "object:map", "object:selection"):
            return False
    return True


def compose(manager, drawings, matrices, disabled):
    """Retain the 3x samples, depth ordering, and renderer-group blend formula."""
    from .export import blend_images, downsample_png
    from .raster import colors, layers
    from .sampling import clipping

    cmd = manager.cmd
    width, height = map(int, matrices[2][2:])
    budget = min((d.sample_budget for d in drawings), default=2048 * 1024**2)
    # Reject oversized backgrounds before PyMOL or NumPy allocates the 3x grid.
    if width * height * 9 * 32 > budget:
        raise ValueError(
            "Ray sample buffers exceed cache_mb; lower dimensions or increase cache_mb"
        )
    enabled = set(cmd.get_names("objects", enabled_only=1))
    for entry in manager.entries.values():
        for name in entry.generated:
            if name in enabled:
                disabled.append(name)
                cmd.disable(name)
    # Let PyMOL supply its background, gradient, and transparent PNG settings.
    cmd.ray(width * 3, height * 3, antialias=0, quiet=1)
    with Image.open(BytesIO(cmd.png(None, prior=1, quiet=1))) as image:
        base = np.array(image.convert("RGBA"))
    nearest = np.full(base.shape[:2], np.inf, np.float32)
    near, far = clipping(cmd.get_view())

    def sample(drawing, meshes):
        rgba, depth = colors(drawing, meshes, matrices, drawing.sample_budget)
        mask = (rgba[:, :, 3] > 0) & (depth >= near) & (depth <= far)
        pixels = np.uint8(np.clip(np.floor(rgba * 255 + 0.5), 0, 255))
        return pixels, np.where(mask, depth, np.inf)

    for drawing in drawings:
        for opacity, meshes in layers(drawing).items():
            if opacity < 0.999999:
                continue
            pixels, depth = sample(drawing, meshes)
            visible = depth < nearest
            base[visible] = pixels[visible]
            nearest[visible] = depth[visible]

    def encode(pixels):
        stream = BytesIO()
        Image.fromarray(pixels).save(stream, format="PNG")
        return downsample_png(stream.getvalue(), width, height)

    result = encode(base)
    passes = []
    for drawing in drawings:
        for opacity, meshes in layers(drawing).items():
            if opacity >= 0.999999:
                continue
            pixels, depth = sample(drawing, meshes)
            visible = depth < nearest
            layer = base.copy()
            layer[visible] = pixels[visible]
            passes.append((opacity, encode(layer)))
    return blend_images(result, passes) if passes else result
