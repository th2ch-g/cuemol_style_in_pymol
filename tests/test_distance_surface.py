"""Reference-scale SES geometry, ownership, and native allocation checks."""

import numpy as np
import pytest
from scipy.spatial import cKDTree

from cuemol_style_in_pymol._edtsurf import distance_surface


def test_distance_surface_preserves_single_atom_radius_and_outward_normals():
    vertices, normals, faces, owners = distance_surface([[0, 0, 0]], [1], detail=10)
    radius = np.linalg.norm(vertices, axis=1)
    assert abs(radius.mean() - 1.7) < 0.1
    assert np.max(abs(radius - 1.7)) < 0.2
    assert np.mean(np.sum(vertices * normals, axis=1) > 0) > 0.99
    np.testing.assert_allclose(np.linalg.norm(normals, axis=1), 1, atol=1e-6)
    assert faces.min() >= 0 and faces.max() < len(vertices)
    assert set(owners) == {0}


def test_disconnected_surfaces_keep_atom_ownership_and_translation():
    positions = np.array([[0, 0, 0], [10, 0, 0]], dtype=np.float32)
    v, _, f, owners = distance_surface(positions, [1, 5])
    assert set(owners) == {0, 1}
    np.testing.assert_array_equal(
        owners, np.argmin(np.linalg.norm(v[:, None] - positions, axis=2), axis=1)
    )
    assert np.linalg.norm(np.ptp(v[f], axis=1), axis=1).max() < 2
    moved = distance_surface(positions + [32, -16, 8], [1, 5])[0] - [32, -16, 8]
    assert cKDTree(v).query(moved)[0].max() < 1e-5
    for left, right in zip(
        distance_surface(positions, [1, 5]), distance_surface(positions, [1, 5])
    ):
        np.testing.assert_array_equal(left, right)


@pytest.mark.parametrize(
    "position, elements, kwargs",
    [
        ([[np.nan, 0, 0]], [1], {}),
        ([[0, 0, 0]], [7], {}),
        ([[0, 0, 0]], [1], {"detail": 0}),
        ([[0, 0, 0]], [1], {"probe": -1}),
        ([[0, 0, 0]], [1], {"max_bytes": 1}),
        ([[0, 0, 0]], [], {}),
    ],
)
def test_invalid_surface_arguments_fail_before_native_allocation(
    position, elements, kwargs
):
    with pytest.raises(ValueError):
        distance_surface(position, elements, **kwargs)
