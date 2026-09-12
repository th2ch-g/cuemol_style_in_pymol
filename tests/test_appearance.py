"""Reference numerical samples and camera-dependent rendering invariants."""

import numpy as np
import pytest
from cuemol_style_in_pymol._edtsurf import pencil, surface

from cuemol_style_in_pymol.geometry import sphere
from cuemol_style_in_pymol.materials import bake, native_blend
from cuemol_style_in_pymol.mesh import Mesh
from cuemol_style_in_pymol.sampling import (
    contour_bands,
    contour_colors,
    material_samples,
    project,
    subdivide,
)


def test_pencil_samples_match_pinned_umbreon_math():
    # Independently evaluated with Umbreon bf75c8a hatch_ink.hpp at 3x.
    samples = np.array(
        [
            [0, 0, 0.2],
            [12.25, 5.75, 0.56],
            [211.5, 37.5, 0.91],
            [-7.1, 140.7, 0.1],
            [140, 300, 0.7],
        ]
    )
    expected = [
        [0.602169514, 0.622299552, 0.611168861],
        [0.606477439, 0.651897013, 0.664469004],
        [0.922826588, 0.913035452, 0.860237598],
        [0.280692011, 0.315371275, 0.334691703],
        [0.464419305, 0.577567399, 0.654062510],
    ]
    actual = pencil(samples[:, :2], samples[:, 2], np.tile([0.25, 0.5, 0.75], (5, 1)))
    np.testing.assert_allclose(actual, expected, atol=2e-7)


def test_pov_toon_and_metal_finishes_are_distinct():
    mesh = sphere([0, 0, 0], 1, [0.8, 0.4, 0.2], 0, 32)
    for a, b in (("toon1", "toon2"), ("diff_metal", "spec_metal")):
        assert np.abs(bake(mesh, a) - bake(mesh, b)).max() > 0.2
    black = bake(mesh, "spec_metal", background=[0, 0, 0])
    white = bake(mesh, "spec_metal", background=[1, 1, 1])
    np.testing.assert_allclose(white, np.clip(black + 0.65, 0, 1), atol=1e-6)


def test_surface_translation_and_voxel_atom_owners():
    coordinates = np.array([[0, 0, 0], [2.3, 0.2, 0], [1.1, 2.1, 0.8]])
    v, n, faces, owners = surface(coordinates, [1, 3, 2])
    moved = surface(coordinates + [16, -8, 4], [1, 3, 2])
    replay = surface(coordinates, [1, 3, 2])
    for actual, expected in zip(replay, (v, n, faces, owners)):
        np.testing.assert_array_equal(actual, expected)
    from scipy.spatial import cKDTree

    # Voxel ties can change under float32 translation by less than half a cell.
    translated = moved[0] - [16, -8, 4]
    assert cKDTree(v).query(translated)[0].max() < 0.25
    assert cKDTree(translated).query(v)[0].max() < 0.25
    assert set(owners) == {0, 1, 2}
    assert len(faces) > 100
    assert np.isfinite(n).all()
    phosphorus = surface([[0, 0, 0]], [5])
    assert len(phosphorus[0]) > 0
    with pytest.raises(RuntimeError, match="budget"):
        surface(coordinates, [1, 3, 2], max_bytes=1)
    with pytest.raises(ValueError, match="finite"):
        surface([[float("nan"), 0, 0]], [1])


def test_native_samples_obey_pixel_bound_and_preserve_interpolation():
    vertices = np.array([[-0.2, -0.2, 0], [0.2, -0.2, 0], [-0.2, 0.2, 0]])
    mesh = Mesh(vertices, [[0, 0, 1]] * 3, vertices + 0.5, [[0, 1, 2]], [0, 1, 2], 0.6)
    matrices = np.eye(4), np.eye(4), np.array([0, 0, 32, 32])
    sampled = subdivide(mesh, matrices, 16 * 1024**2)
    pixels = project(sampled.vertices, matrices)[0][sampled.faces]
    assert (
        np.linalg.norm(pixels - np.roll(pixels, -1, axis=1), axis=2).max()
        <= 1 / 3 + 1e-6
    )
    np.testing.assert_allclose(sampled.colors, sampled.vertices + 0.5, atol=1e-6)
    assert sampled.opacity == 0.6
    with pytest.raises(ValueError, match="cache_mb"):
        subdivide(mesh, matrices, 16)


def test_sampling_clips_crossing_triangles_before_subdivision():
    mesh = Mesh(
        [[-2, 0, 0], [0, -2, 0], [0, 0, 2]],
        [[0, 0, 1]] * 3,
        [[1, 0, 0]] * 3,
        [[0, 1, 2]],
        [0, 1, 2],
    )
    matrices = np.eye(4), np.eye(4), np.array([0, 0, 8, 8])
    result = subdivide(mesh, matrices, 16 * 1024**2)
    assert len(result.faces)
    assert np.abs(result.vertices).max() <= 1 + 1e-6


def test_transparent_group_samples_only_the_front_surface():
    front = np.array([[-0.4, -0.4, -0.2], [0.4, -0.4, -0.2], [0, 0.4, -0.2]])
    vertices = np.vstack((front, front + [0, 0, 0.4]))
    mesh = Mesh(
        vertices,
        [[0, 0, 1]] * 6,
        [[0.3, 0.6, 0.8]] * 6,
        [[0, 1, 2], [3, 4, 5]],
        [0, 0, 0, 1, 1, 1],
        0.6,
    )
    matrices = np.eye(4), np.eye(4), np.array([0, 0, 48, 48])
    visible, colors = material_samples(
        mesh, "default", matrices, (1, 1, 1), 16 * 1024**2
    )
    assert len(visible.faces) > 100
    assert set(visible.owners[visible.faces].ravel()) == {0}
    assert np.isfinite(colors).all()
    assert mesh.opacity == 0.6
    assert 0.6 < visible.opacity < 1


def test_native_group_alpha_matches_reference_display_blend():
    # CueMol's 75%-opaque surface sample is (255, 34, 34), not (255, 70, 70).
    alpha, colors = native_blend([[1, 9 / 255, 9 / 255]], 0.75, (1, 1, 1))
    composite = (alpha * colors + 1 - alpha) * 255
    np.testing.assert_allclose(composite, [[255, 34, 34]], atol=1)
    for background in ((0, 0, 0), (1, 1, 1), (0.2, 0.7, 0.5)):
        for opacity in (0, 0.15, 0.75, 1):
            alpha, colors = native_blend([[0, 0, 0], [1, 1, 1]], opacity, background)
            assert 0 <= alpha <= 1
            assert ((colors >= 0) & (colors <= 1)).all()
            if opacity == 0:
                np.testing.assert_allclose(
                    alpha * colors + (1 - alpha) * np.array(background),
                    [background, background],
                    atol=1e-7,
                )


def test_native_contour_samples_keep_ink_on_the_front_surface():
    mesh = Mesh(
        [[-0.5, 0, -2], [0.5, 0, -2], [0, 0.5, -2]],
        [[0, 0, 1]] * 3,
        [[0.8, 0.8, 0.8]] * 3,
        [[0, 1, 2]],
        [0, 0, 0],
    )
    projection = np.eye(4)
    projection[2, 2] = -0.1
    matrices = np.eye(4), projection, np.array([0, 0, 48, 48])
    edges = np.array([[-0.5, 0, -2, 0.5, 0, -2, 0, 0, 1, 0, 0, -1]])
    colors = contour_colors(
        mesh.colors.copy(), mesh, matrices, edges, 0.06, (0, 0, 0), (1, 1, 1), (2, 10)
    )
    np.testing.assert_array_equal(colors[:2], 0)
    np.testing.assert_array_equal(colors[2], mesh.colors[2])


def test_png_sample_average_preserves_alpha_and_ignores_hidden_rgb():
    from io import BytesIO

    from PIL import Image

    from cuemol_style_in_pymol.export import downsample_png

    pixels = np.full((3, 6, 4), [0, 0, 255, 0], dtype=np.uint8)
    pixels[1, 1] = [255, 0, 0, 255]
    pixels[:, 3:] = [20, 80, 160, 255]
    output = BytesIO()
    Image.fromarray(pixels).save(output, "PNG")
    result = Image.open(BytesIO(downsample_png(output.getvalue(), 2, 1)))
    np.testing.assert_array_equal(
        np.asarray(result), [[[255, 0, 0, 28], [20, 80, 160, 255]]]
    )
    with pytest.raises(RuntimeError, match="size"):
        downsample_png(output.getvalue(), 3, 1)


@pytest.mark.parametrize("perspective", [False, True])
def test_ray_contours_keep_screen_width_and_do_not_show_hidden_edges(perspective):
    mesh = Mesh(
        [[-0.9, -0.9, -2], [0.9, -0.9, -2], [0.9, 0.9, -2], [-0.9, 0.9, -2]],
        [[0, 0, 1]] * 4,
        [[1, 0, 0]] * 4,
        [[0, 1, 2], [0, 2, 3]],
        [0] * 4,
    )
    projection = np.eye(4)
    projection[2, 2] = -0.1
    if perspective:
        projection[0, 0] = projection[1, 1] = 2
        projection[2, 2], projection[2, 3] = -11 / 9, -20 / 9
        projection[3, 2], projection[3, 3] = -1, 0
    matrices = np.eye(4), projection, np.array([0, 0, 48, 48])
    edges = np.array([[-0.5, 0, -2, 0.5, 0, -2, 0, 0, 1, 0, 0, -1]], float)
    args = matrices, 0.02, (0, 0, 0), (1, 1, 1), (2, 10), 16 * 1024**2
    bands = contour_bands(mesh, edges, *args)
    assert len(bands.faces)
    pixels = project(bands.vertices, matrices)[0]
    assert np.ptp(pixels[:, 1]) == pytest.approx(1, abs=1e-4)
    edges[:, [2, 5]] = -3
    assert not len(contour_bands(mesh, edges, *args).faces)
