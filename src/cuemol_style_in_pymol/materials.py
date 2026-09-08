"""Material colors for native transparency and ray export."""

import numpy as np

from .mesh import unit
from .presets import MATERIALS


def bake(mesh, material, rotation=None):
    """Approximate the GLSL material using view-space vertex samples."""
    rotation = np.eye(3) if rotation is None else np.asarray(rotation).reshape(3, 3)
    n = unit(mesh.normals @ rotation.T)
    light = unit([0.35, 0.65, 1.0])
    diffuse = np.maximum(n @ light, 0)
    specular = np.maximum(n @ unit(light + [0, 0, 1]), 0)
    color = mesh.colors.copy()
    if material == "nolighting":
        return color
    shade = 0.2 + 0.8 * diffuse
    if material == "toon1":
        shade = np.where(diffuse < 0.2, 0.48, np.where(diffuse < 0.65, 0.76, 1.0))
    elif material == "toon2":
        shade = np.where(diffuse < 0.5, 0.5, 1.0)
    elif material == "shadow":
        shade = np.full(len(n), 0.75)
    elif material == "matte":
        shade = 0.3 + 0.6 * diffuse
    elif material in ("diff_metal", "spec_metal"):
        shade = 0.2 + 0.5 * diffuse
    elif material in ("metallic_chrome", "metallic_copper"):
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
    if material in (
        "diff_metal",
        "spec_metal",
        "metallic_chrome",
        "metallic_copper",
    ):
        power = 70.0 if material in ("spec_metal", "metallic_chrome") else 25.0
        result += 0.7 * specular[:, None] ** power
    return np.clip(result, 0, 1)


def material_id(name):
    return MATERIALS.index(name)
