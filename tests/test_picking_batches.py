"""Ray culling must retain thin surfaces, partial batches, and exact owners."""

import numpy as np
import pytest

from cuemol_style_in_pymol.mesh import Mesh, ray_hits


def triangles():
    offsets = np.c_[np.arange(130) * 3, np.zeros((130, 2))]
    vertices = (
        offsets[:, None] + [[-0.4, -0.4, 0], [0.4, -0.4, 0], [0, 0.4, 0]]
    ).reshape(-1, 3)
    return Mesh(
        vertices,
        np.tile([0, 0, 1], (390, 1)),
        np.ones((390, 3)),
        np.arange(390).reshape(-1, 3),
        np.repeat(np.arange(130), 3),
    )


@pytest.mark.parametrize("face", [0, 63, 64, 65, 127, 129])
def test_axis_parallel_ray_keeps_partial_batches_and_owners(face):
    mesh = triangles()
    hit = ray_hits(mesh, np.array([face * 3, 0, 2]), np.array([0, 0, -1]))
    assert hit == (2.0, face)


def test_misses_and_away_facing_rays_do_not_pick():
    mesh = triangles()
    assert ray_hits(mesh, np.array([1, 0, 2]), np.array([0, 0, -1])) is None
    assert ray_hits(mesh, np.array([0, 0, -2]), np.array([0, 0, -1])) is None
    assert ray_hits(mesh, np.array([0, 0, 2]), np.array([1, 0, 0])) is None


def test_oblique_ray_and_empty_geometry():
    assert ray_hits(triangles(), np.array([-2, 0, 2]), np.array([1, 0, -1])) == (2.0, 0)
    assert ray_hits(Mesh([], [], [], [], []), np.zeros(3), np.ones(3)) is None
