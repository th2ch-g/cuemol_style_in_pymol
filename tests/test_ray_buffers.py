"""Reject oversized ray grids before allocating native or NumPy image buffers."""

from types import SimpleNamespace

import numpy as np
import pytest

from cuemol_style_in_pymol.ray import compose


def test_oversized_image_fails_before_entering_pymol():
    manager = SimpleNamespace(cmd=object())
    drawings = [SimpleNamespace(sample_budget=1024 * 1024)]
    matrices = np.eye(4), np.eye(4), np.array([0, 0, 32768, 32768])
    with pytest.raises(ValueError, match="Ray sample buffers exceed cache_mb"):
        compose(manager, drawings, matrices, [])
