"""Material colors for native transparency and ray export."""

import numpy as np

from .mesh import unit
from .presets import MATERIALS, OPENGL_MATERIALS


def drawing_tone(normals):
    """Camera-space tone, separate from the pigment used for pencil marks."""
    n = unit(normals)
    n = n * np.where(n[:, 2:3] < 0, -1, 1)
    flash = np.where(n[:, 2] > 0, np.clip((n[:, 2] + 0.5) / 1.5, 0, 1), 0)
    ndl = n @ unit([1, 1, 1])
    key = np.where(ndl > 0, np.clip((ndl + 0.5) / 1.5, 0, 1), 0)
    tone = 0.05 + 0.85 * 1.302 * (0.6 * flash + 0.4 * key)
    tone *= 1 - (1 - np.maximum(n[:, 2], 0)) ** 3.5 * (1 - 0.35 * np.clip(tone, 0, 1))
    tone = np.clip(tone / 1.2, 0, 1) ** 2.4
    tone = np.where(tone <= 0.0031308, 12.92 * tone, 1.055 * tone ** (1 / 2.4) - 0.055)
    knee = np.clip((tone - 0.81) / 0.05, 0, 1)
    return tone + (1 - tone) * knee**2 * (3 - 2 * knee)


def pencil_average(colors, tone):
    """Filtered pencil coverage for native vertex-color-only rendering.

    Native CGO cannot evaluate fragment shaders. Averaging the marks retains
    the paper and pigment tone without an aliased, camera-frozen pixel grid.
    """
    paper = np.array([0.941, 0.925, 0.867])
    ink = colors * (0.4 + 0.6 * tone[:, None])
    luma = np.array([0.2126, 0.7152, 0.0722])
    ink *= np.minimum(1, (paper @ luma) * 0.85 / np.maximum(ink @ luma, 1e-8))[:, None]
    result = np.tile(paper, (len(colors), 1))
    for threshold, strength in ((0.92, 1.0), (0.62, 0.74), (0.34, 0.38)):
        growth = np.clip((threshold - tone) * 10, 0, 1)
        # Reference 3x ink integration: full-layer mean coverage is 0.256.
        coverage = growth * (0.275 - 0.019 * growth)
        result *= 1 - coverage[:, None] * (1 - ink * strength)
    return np.clip(result, 0, 1)


def bake(mesh, material, rotation=None):
    """Approximate the GLSL material using view-space vertex samples."""
    rotation = np.eye(3) if rotation is None else np.asarray(rotation).reshape(3, 3)
    n = unit(mesh.normals @ rotation.T)
    light = unit([1.0, 1.0, 1.5])
    diffuse = np.maximum(n @ light, 0)
    specular = np.maximum(n @ unit(light + [0, 0, 1]), 0)
    color = mesh.colors.copy()
    if material == "richardson":
        return pencil_average(color, drawing_tone(n))
    ambient, lambert, spec, power = material_coefficients(material)
    shade = ambient + lambert * diffuse
    if material in ("metallic_chrome", "metallic_copper"):
        bands = 0.3 + 0.7 * (0.5 + 0.5 * np.sin(n[:, 1] * 11.0 + n[:, 2] * 4.0)) ** 2
        metal = (
            [1.0, 0.55, 0.27] if material == "metallic_copper" else [0.83, 0.9, 0.96]
        )
        color = np.tile(metal, (len(n), 1)) * (0.65 + 0.35 * color)
        shade = bands
    elif material == "stone35":
        p = mesh.vertices
        grain = np.sin(p @ [12.1, 17.3, 8.7]) * np.sin(p @ [29.7, 4.2, 21.3])
        color *= 0.65 + 0.35 * grain[:, None]
    elif material in ("wood31", "wood14scl2"):
        p = mesh.vertices
        ring = np.linalg.norm(p[:, [0, 2]], axis=1)
        grain = 0.5 + 0.5 * np.sin(
            ring * (3.0 if material == "wood31" else 7.0) + 0.5 * np.sin(p[:, 1])
        )
        color = np.array([0.25, 0.095, 0.035]) + grain[:, None] * [0.55, 0.35, 0.15]
    result = color * shade[:, None]
    if material in ("metallic_chrome", "metallic_copper"):
        spec, power = 0.7, 70.0 if material == "metallic_chrome" else 25.0
    result += spec * np.where(diffuse > 0, specular**power, 0)[:, None]
    return np.clip(result, 0, 1)


def material_id(name):
    return MATERIALS.index(name)


def material_coefficients(name):
    return OPENGL_MATERIALS.get(name, OPENGL_MATERIALS["default"])
