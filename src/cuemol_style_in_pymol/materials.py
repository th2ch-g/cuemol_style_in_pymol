"""Material colors for native transparency and ray export."""

import numpy as np

from .mesh import unit
from .presets import MATERIALS, PBR_MATERIALS, POV_MATERIALS


def srgb_encode(values):
    values = np.clip(values, 0, 1)
    return np.where(
        values <= 0.0031308, 12.92 * values, 1.055 * values ** (1 / 2.4) - 0.055
    )


def srgb_decode(values):
    values = np.maximum(values, 0)
    return np.where(
        values <= 0.04045, values / 12.92, ((values + 0.055) / 1.055) ** 2.4
    )


def native_opacity(opacity, background):
    """Choose a native alpha that can represent every encoded blend color."""
    if opacity == 0 or opacity == 1:
        return float(opacity)
    background = np.asarray(background)
    encoded = srgb_encode(background)
    dark = 1 - srgb_decode((1 - opacity) * encoded) / np.maximum(background, 1e-12)
    light = (srgb_decode(opacity + (1 - opacity) * encoded) - background) / np.maximum(
        1 - background, 1e-12
    )
    return float(
        np.clip(
            max(
                opacity,
                np.max(np.where(background > 0, dark, 0)),
                np.max(np.where(background < 1, light, 0)),
            ),
            0,
            1,
        )
    )


def native_blend(colors, opacity, background):
    """Match Umbreon's sRGB group blend over a solid background using CGO."""
    background = np.asarray(background)
    alpha = native_opacity(opacity, background)
    blended = srgb_decode(
        opacity * srgb_encode(colors) + (1 - opacity) * srgb_encode(background)
    )
    return alpha, np.clip(
        (blended - (1 - alpha) * background) / max(alpha, 1e-12), 0, 1
    )


def drawing_tone(normals, view_direction=None, fog=None):
    """Camera-space tone, separate from the pigment used for pencil marks."""
    n = unit(normals)
    view = (
        np.tile([0, 0, 1], (len(n), 1))
        if view_direction is None
        else unit(view_direction)
    )
    n = n * np.where(np.sum(n * view, axis=1, keepdims=True) < 0, -1, 1)
    flash = np.where(n[:, 2] > 0, np.clip((n[:, 2] + 0.5) / 1.5, 0, 1), 0)
    ndl = n @ unit([1, 1, 1])
    key = np.where(ndl > 0, np.clip((ndl + 0.5) / 1.5, 0, 1), 0)
    tone = 0.05 + 0.85 * 1.302 * (0.6 * flash + 0.4 * key)
    tone *= 1 - (1 - np.maximum(np.sum(n * view, axis=1), 0)) ** 3.5 * (
        1 - 0.35 * np.clip(tone, 0, 1)
    )
    if fog is not None:
        tone = 1 - (1 - tone) * fog
    tone = np.clip(tone / 1.2, 0, 1) ** 2.4
    tone = np.where(tone <= 0.0031308, 12.92 * tone, 1.055 * tone ** (1 / 2.4) - 0.055)
    knee = np.clip((tone - 0.81) / 0.05, 0, 1)
    return tone + (1 - tone) * knee**2 * (3 - 2 * knee)


def pencil_average(colors, tone):
    """Filtered pencil coverage for native vertex-color-only rendering.

    Native CGO cannot evaluate fragment shaders. Averaging the marks retains
    the paper and pigment tone without an aliased, camera-frozen pixel grid.
    """
    paper = np.array([240, 236, 221]) / 255
    ink = colors * (0.4 + 0.6 * tone[:, None])
    luma = np.array([0.2126, 0.7152, 0.0722])
    ink *= np.minimum(1, (paper @ luma - 0.15) / np.maximum(ink @ luma, 1e-8))[:, None]
    result = np.tile(paper, (len(colors), 1))
    for threshold, strength in ((0.92, 1.0), (0.62, 0.74), (0.34, 0.38)):
        growth = np.clip((threshold - tone) * 10, 0, 1)
        # Reference 3x ink integration: full-layer mean coverage is 0.256.
        coverage = growth * (0.275 - 0.019 * growth)
        result *= 1 - coverage[:, None] * (1 - ink * strength)
    return np.clip(result, 0, 1)


def bake(
    mesh, material, rotation=None, background=(0, 0, 0), view_direction=None, fog=None
):
    """Evaluate the current CueMol direct material at native CGO samples."""
    rotation = np.eye(3) if rotation is None else np.asarray(rotation).reshape(3, 3)
    n = unit(mesh.normals @ rotation.T)
    view = (
        np.tile([0, 0, 1], (len(n), 1))
        if view_direction is None
        else unit(view_direction)
    )
    n *= np.where(np.sum(n * view, axis=1, keepdims=True) < 0, -1, 1)
    light = unit([1.0, 1.0, 1.0])
    diffuse = np.maximum(n @ light, 0)
    specular = np.maximum(np.sum(n * unit(light + view), axis=1), 0)
    color = mesh.colors.copy()
    if material == "richardson":
        return pencil_average(color, drawing_tone(n))
    if material in PBR_MATERIALS:
        ambient, lambert, metal, roughness, spec, reflection = PBR_MATERIALS[material]
        alpha2 = max(roughness**2, 0.001) ** 2
        ndv = np.maximum(np.sum(n * view, axis=1), 1e-4)
        f0 = 0.08 * spec * (1 - metal) + color * metal
        half = unit(light + view)
        vdoth = np.clip(np.sum(view * half, axis=1), 0, 1)
        fresnel = f0 + (1 - f0) * (1 - vdoth[:, None]) ** 5

        def masking(cosine):
            squared = np.maximum(cosine**2, 1e-12)
            return 0.5 * (np.sqrt(1 + alpha2 * (1 - squared) / squared) - 1)

        distribution = alpha2 / (specular**2 * (alpha2 - 1) + 1) ** 2
        geometry = 1 / (1 + masking(ndv) + masking(diffuse))
        highlight = np.where(diffuse > 0, 0.52 * distribution * geometry / (4 * ndv), 0)
        result = (
            color
            * (
                ambient
                + lambert
                * (1 - metal)
                * (0.52 * diffuse + 0.78 * np.maximum(n[:, 2], 0))
            )[:, None]
        )
        result += highlight[:, None] * fresnel * (f0.max(axis=1) > 0)[:, None]
        result += (reflection if reflection > 0 else f0) * np.asarray(background)
        return np.clip(
            result
            if fog is None
            else np.asarray(background) + (result - background) * fog[:, None],
            0,
            1,
        )
    ambient, lambert, spec, power, brilliance, phong, phong_size, reflection, metal = (
        finish(material)
    )
    flash = np.maximum(n[:, 2], 0)
    shade = ambient + lambert * (
        0.52 * np.where(diffuse > 0, diffuse**brilliance, 0)
        + 0.78 * np.where(flash > 0, flash**brilliance, 0)
    )
    result = color * shade[:, None]
    reflected = 2 * diffuse[:, None] * n - light
    rv = np.maximum(np.sum(reflected * view, axis=1), 0)
    highlight = 0.52 * np.where(
        diffuse > 0, spec * specular**power + phong * rv**phong_size, 0
    )
    tint = np.ones_like(color)
    if metal:
        angle = np.arccos(np.clip(diffuse, 0, 1)) * (2 / np.pi)
        fresnel = np.clip(0.014567225 / (angle - 1.12) ** 2 - 0.011612903, 0, 1)
        tint = fresnel[:, None] + (1 - fresnel[:, None]) * color
    result += highlight[:, None] * tint + reflection * np.asarray(background)
    return np.clip(
        result
        if fog is None
        else np.asarray(background) + (result - background) * fog[:, None],
        0,
        1,
    )


def material_id(name):
    return MATERIALS.index(name)


def material_coefficients(name):
    return finish(name)[:4]


def finish(name):
    return POV_MATERIALS.get(name, POV_MATERIALS["default"])
