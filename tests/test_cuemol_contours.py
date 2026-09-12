"""Regression checks for screen contours and renderer-level transparency."""

from io import BytesIO
from types import SimpleNamespace

import numpy as np
from PIL import Image
from scipy.ndimage import binary_erosion, distance_transform_edt
from scipy.spatial.transform import Rotation

from cuemol_style_in_pymol.contours import contour_image, stroke_image
from cuemol_style_in_pymol.export import blend_images
from cuemol_style_in_pymol.materials import srgb_decode, srgb_encode
from cuemol_style_in_pymol.mesh import Mesh
from cuemol_style_in_pymol.presets import Profile
from cuemol_style_in_pymol.raster import colors, pixel_mesh


def sheet(diagonal):
    x, y = np.meshgrid(np.linspace(-1, 1, 7), np.linspace(-1, 1, 7))
    vertices = np.c_[x.ravel(), y.ravel(), 0.012 * (x * y).ravel()]
    quads = np.array(
        [
            [a, a + 1, a + 8, a + 7]
            for row in range(6)
            for a in range(row * 7, row * 7 + 6)
        ]
    )
    faces = quads[
        :, [[0, 1, 2], [0, 2, 3]] if diagonal else [[0, 1, 3], [1, 2, 3]]
    ].reshape(-1, 3)
    return Mesh(
        vertices,
        [[0, 0, 1]] * len(vertices),
        [[0.2, 0.5, 0.7]] * len(vertices),
        faces,
        np.zeros(len(vertices), int),
    )


def camera(angle):
    view = np.eye(4)
    view[:3, :3] = Rotation.from_euler("y", angle, degrees=True).as_matrix()
    view[2, 3] = -5
    projection = np.diag([0.65, 0.65, -0.1, 1])
    return view, projection, np.array([0, 0, 80, 80])


def test_rotating_sheet_does_not_outline_internal_triangles():
    for angle in (0, 40, 75, 85, 88, 89.5, 90, 90.5, 92, 105, 180):
        for diagonal in (0, 1):
            surface, _, chains, width = contour_image(
                [sheet(diagonal)], camera(angle), 0.03, 2**26
            )
            ink = stroke_image(chains, width, surface.shape[:2]) > 0
            body = surface[:, :, 0] > 0
            boundary = body & ~binary_erosion(body)
            distance = distance_transform_edt(~boundary)
            assert ink.any()
            assert distance[ink].max() <= width + 2
            assert not np.any(ink & (distance_transform_edt(body) > width + 2))


def test_transparent_group_keeps_only_nearest_surface_and_one_ink_layer():
    front, rear = sheet(0), sheet(1)
    rear.vertices[:, 2] -= 1
    rear.colors[:] = [1, 0, 0]
    drawing = SimpleNamespace(
        profile=Profile(material="nolighting", edges="edges"),
        edge_color=(0, 0, 0),
        background=(1, 1, 1),
        fog=(100, 200),
    )
    expected, _ = colors(drawing, [front], camera(0), 2**27)
    actual, depth = colors(drawing, [rear, front], camera(0), 2**27)
    np.testing.assert_array_equal(actual, expected)
    sampled = pixel_mesh(actual, depth, camera(0), 2**27, 0.4, drawing.background)
    assert len(sampled.faces) == 2 * np.count_nonzero(actual[:, :, 3])
    assert np.all(sampled.vertices[:, 2] > -0.02)


def test_completed_group_blend_is_affine_and_order_independent():
    def png(rgb, alpha=255):
        stream = BytesIO()
        Image.new("RGBA", (2, 2), (*rgb, alpha)).save(stream, format="PNG")
        return stream.getvalue()

    base = png((230, 240, 250), 0)
    a, b = png((30, 60, 90)), png((110, 150, 190))
    actual = blend_images(base, [(0.7, a), (0.5, b)])
    assert actual == blend_images(base, [(0.5, b), (0.7, a)])
    values = np.asarray(Image.open(BytesIO(actual)))
    expected = srgb_decode(
        -0.2 * srgb_encode(np.array([230, 240, 250]) / 255)
        + 0.7 * srgb_encode(np.array([30, 60, 90]) / 255)
        + 0.5 * srgb_encode(np.array([110, 150, 190]) / 255)
    )
    np.testing.assert_allclose(
        values[0, 0, :3], np.clip(np.floor(expected * 255 + 0.5), 0, 255), atol=1
    )
    assert values[0, 0, 3] == 255
