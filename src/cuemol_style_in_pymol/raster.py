"""Shared visible-surface samples for contours and native CGO output."""

import numpy as np

from ._edtsurf import pencil
from .contours import contour_image, stroke_image, surface_image
from .materials import bake, drawing_tone, native_blend
from .mesh import Mesh, merge


def layers(drawing):
    return {
        opacity: [p.mesh for p in drawing.pieces if p.mesh.opacity == opacity]
        for opacity in sorted({p.mesh.opacity for p in drawing.pieces})
        if opacity > 0
    }


def buffers(drawing, meshes, matrices, budget):
    if drawing.profile.edges == "none":
        surface, depth = surface_image(merge(meshes), matrices, budget)
        return surface, depth, np.zeros(depth.shape, np.float32)
    surface, depth, chains, width = contour_image(
        meshes,
        matrices,
        drawing.profile.edge_width,
        budget,
        drawing.profile.edges == "silhouette",
        drawing.fog[0] + (drawing.fog[1] - drawing.fog[0]) * 0.2,
        drawing.profile.representation in ("cpk", "ballstick", "sticks", "nucleic"),
        drawing.fog[0],
    )
    ink = stroke_image(chains, width, depth.shape)
    # The reference permits a one-output-pixel halo at an occluding contour.
    if width > 6:
        from scipy.ndimage import distance_transform_edt

        backbone = stroke_image(chains, 1, depth.shape) > 0
        beyond_halo = distance_transform_edt(~backbone) > 3
        pixel = 2 / (matrices[1][1, 1] * depth.shape[0])
        if abs(matrices[1][3, 3]) < 0.5:
            pixel = pixel * np.maximum(surface[:, :, 0], 1e-8)
        hidden = (
            beyond_halo & (surface[:, :, 0] > 0) & (ink - surface[:, :, 0] > 12 * pixel)
        )
        ink[hidden] = 0
    return surface, depth, ink


def colors(drawing, meshes, matrices, budget):
    surface, _, ink = buffers(drawing, meshes, matrices, budget)
    mask = surface[:, :, 0] > 0
    h, w = mask.shape
    yy, xx = np.nonzero(mask)
    direction = None
    projection = matrices[1]
    if abs(projection[3, 3]) < 0.5:
        direction = np.c_[
            -((2 * (xx + 0.5) / w - 1) + projection[0, 2]) / projection[0, 0],
            -((1 - 2 * (yy + 0.5) / h) + projection[1, 2]) / projection[1, 1],
            np.ones(len(xx)),
        ]
    fog = np.clip(
        (drawing.fog[1] - surface[mask, 0])
        / max(drawing.fog[1] - drawing.fog[0], 1e-8),
        0,
        1,
    )
    result = np.zeros((h, w, 4), np.float32)
    ink_mask = ink > 0
    edge_fog = np.clip(
        (drawing.fog[1] - ink[ink_mask]) / max(drawing.fog[1] - drawing.fog[0], 1e-8),
        0,
        1,
    )
    background = np.asarray(drawing.background)
    ink_color = (
        background + (np.asarray(drawing.edge_color) - background) * edge_fog[:, None]
    )
    if drawing.profile.material == "richardson":
        tone = drawing_tone(surface[mask, 1:4], direction, fog)
        base = np.empty((h, w, 3), np.float32)
        base[:] = np.array([240, 236, 221]) / 255
        base[ink_mask] = ink_color
        result[mask, :3] = pencil(
            np.c_[(xx + 0.5) / 3, (yy + 0.5) / 3],
            tone,
            surface[mask, 4:7],
            False,
            base[mask],
        )
        result[ink_mask & ~mask, :3] = base[ink_mask & ~mask]
    else:
        samples = Mesh(
            np.zeros((len(xx), 3)), surface[mask, 1:4], surface[mask, 4:7], [], []
        )
        result[mask, :3] = bake(
            samples,
            drawing.profile.material,
            background=background,
            view_direction=direction,
            fog=fog,
        )
        result[ink_mask, :3] = ink_color
    result[mask | ink_mask, 3] = 1
    view_z = surface[:, :, 0].copy()
    view_z[ink_mask] = np.where(
        mask[ink_mask], np.minimum(view_z[ink_mask], ink[ink_mask]), ink[ink_mask]
    )
    return result, view_z


def pixel_mesh(image, view_z, matrices, budget, opacity=1, background=(0, 0, 0)):
    """One depth-tested quad per covered sample, with no overlapping back faces."""
    yy, xx = np.nonzero(image[:, :, 3] > 0)
    if len(xx) * 400 > budget:
        raise ValueError("Visible pixel CGO exceeds cache_mb")
    h, w = view_z.shape
    rgb = image[yy, xx, :3]
    if opacity < 1:
        opacity, rgb = native_blend(rgb, opacity, background)
    offsets = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    xy = np.c_[xx, yy][:, None, :] + offsets
    projection = matrices[1]
    eye_z = -view_z[yy, xx, None]
    clip_w = projection[3, 2] * eye_z + projection[3, 3]
    eye = np.empty((len(xx), 4, 3))
    eye[:, :, 0] = (
        (2 * xy[:, :, 0] / w - 1) * clip_w - projection[0, 2] * eye_z - projection[0, 3]
    ) / projection[0, 0]
    eye[:, :, 1] = (
        (1 - 2 * xy[:, :, 1] / h) * clip_w - projection[1, 2] * eye_z - projection[1, 3]
    ) / projection[1, 1]
    eye[:, :, 2] = eye_z
    inverse = np.linalg.inv(matrices[0])
    vertices = eye.reshape(-1, 3) @ inverse[:3, :3].T + inverse[:3, 3]
    faces = (np.arange(len(xx))[:, None, None] * 4 + [[0, 2, 1], [0, 3, 2]]).reshape(
        -1, 3
    )
    normals = np.tile(matrices[0][2, :3], (len(vertices), 1))
    return Mesh(
        vertices,
        normals,
        np.repeat(rgb, 4, axis=0),
        faces,
        np.zeros(len(vertices), int),
        opacity,
    )
