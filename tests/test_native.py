"""Compact geometry must retain atom identity, depth, state, and lifecycle."""

from io import BytesIO
from unittest.mock import patch

import numpy as np
import pymol2
import pytest
from PIL import Image

from cuemol_style_in_pymol import cuemol_style
from cuemol_style_in_pymol.controller import manager_for
from cuemol_style_in_pymol.native import Builder
from cuemol_style_in_pymol.source import read


def sample(cmd):
    from chempy import Atom, Bond
    from chempy.models import Indexed

    model = Indexed()
    for index, symbol in enumerate(("C", "C", "O", "N", "H", "H")):
        atom = Atom()
        atom.name, atom.symbol = f"{symbol}{index + 1}", symbol
        atom.resn, atom.resi = "LIG", "1"
        atom.coord = [index * 1.3, index % 2 * 0.3, 0]
        model.atom.append(atom)
        if index:
            bond = Bond()
            bond.index, bond.order = [index - 1, index], 1
            model.bond.append(bond)
    cmd.load_model(model, "sample", zoom=0)


def apply(cmd, **kwargs):
    return cuemol_style(
        "ballstick", quality="low", atomic_mode="native", quiet=1, _self=cmd, **kwargs
    )


def reset(cmd):
    cuemol_style("reset", name="all", _self=cmd)


def test_analytic_depth_and_owner_for_spheres_and_capped_cylinders():
    builder = Builder(("sphere", "cylinder"))
    builder.sphere([0, 0, 0], 1, [1, 0, 0], 0, 8)
    builder.cylinder([3, 0, -1], [3, 0, 1], 0.5, [0, 1, 0], 1, 8)
    native = builder.finish()
    assert native.hit(np.array([0, 0, 3]), np.array([0, 0, -1])) == (2, 0)
    assert native.hit(np.array([0, 0, 0]), np.array([0, 0, 1])) == (1, 0)
    assert native.hit(np.array([3, 0, 3]), np.array([0, 0, -1])) == (2, 1)
    assert native.hit(np.array([5, 0, 0]), np.array([-1, 0, 0])) == (1.5, 1)
    assert native.hit(np.array([5, 0, 3]), np.array([1, 0, 0])) is None
    assert len(native.cgo()) * 4 == native.cgo_nbytes
    native.opacity = 0
    assert native.hit(np.array([0, 0, 3]), np.array([0, 0, -1])) is None


@pytest.mark.parametrize("opacity", [0, 0.4])
def test_native_states_refresh_ray_and_reset(opacity, tmp_path):
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        sample(cmd)
        original_mode = cmd.get_setting_int("transparency_mode")
        cmd.show_as("sticks", "sample")
        coordinates = cmd.get_coords("sample").copy()
        cmd.create("sample", "sample", source_state=1, target_state=2, zoom=0)
        cmd.translate([0, 0, 2], "sample", state=2, camera=0)
        before = cmd.get_model("sample")
        _, saved = read(cmd, "sample")
        entry = apply(cmd, transparency=opacity)
        assert cmd.get_setting_int("transparency_mode") == (
            3 if opacity else original_mode
        )
        drawings = next(iter(entry.drawings.values()))
        assert len(drawings) == 2
        assert len(drawings[0].native) == 1
        part = drawings[0].native[0]
        assert len(part.spheres) == len(before.atom)
        assert part.opacity == pytest.approx(1 - opacity)
        assert part.cgo_nbytes < len(before.atom) * 200
        np.testing.assert_array_equal(part.spheres[:, :3], coordinates)
        np.testing.assert_allclose(
            drawings[1].native[0].spheres[:, :3], coordinates + [0, 0, 2], atol=1e-6
        )
        assert all(cmd.get_type(name) == "object:cgo" for name in entry.native_objects)
        cmd.orient("sample")
        image = Image.open(BytesIO(cmd.png(None, 160, 120, ray=1, quiet=1)))
        assert np.ptp(np.asarray(image)[..., :3]) > 100
        output = tmp_path / "native-ray.png"
        cuemol_style("ray", filename=str(output), width=160, height=120, _self=cmd)
        assert np.ptp(np.asarray(Image.open(output))[..., :3]) > 100
        np.testing.assert_array_equal(cmd.get_coords("sample", state=1), coordinates)
        cmd.translate([1, 0, 0], "sample", state=1, camera=0)
        cuemol_style("refresh", _self=cmd)
        refreshed = next(iter(manager_for(cmd).entries.values()))
        np.testing.assert_allclose(
            next(iter(refreshed.drawings.values()))[0].native[0].spheres[:, :3],
            coordinates + [1, 0, 0],
            atol=1e-6,
        )
        reset(cmd)
        _, restored = read(cmd, "sample")
        assert restored == saved
        assert cmd.get_names("objects") == ["sample"]
        assert cmd.get_setting_int("transparency_mode") == original_mode


def test_subset_isolation_and_palette_queries_scale_with_distinct_colors():
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        sample(cmd)
        with patch.object(
            cmd, "get_color_tuple", wraps=cmd.get_color_tuple
        ) as get_color:
            snapshots, _ = read(cmd, "sample")
        state = snapshots["sample"][0]
        assert get_color.call_count == len({a.color for a in state.atoms})
        mask = np.arange(len(state.atoms)) < 4
        selected = state.subset(mask)
        before = list(state.model.atom[0].coord)
        selected.model.atom[0].coord[0] += 10
        selected.model.atom[0].name = "changed"
        assert state.model.atom[0].coord == before
        assert state.model.atom[0].name != "changed"
        assert all(max(b.index) < 4 for b in selected.model.bond)
        shared = state.subset(mask, copy_atoms=False)
        assert shared.model.atom[0] is state.model.atom[0]
        assert shared.model.bond[0] is not state.model.bond[0]


def test_invalid_native_options_preserve_previous_view():
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        sample(cmd)
        entry = apply(cmd)
        with pytest.raises(Exception, match="atomic_mode"):
            cuemol_style(atomic_mode="invalid", _self=cmd)
        assert manager_for(cmd).entries[entry.name] is entry
        with pytest.raises(Exception, match="Native atoms require"):
            cuemol_style("richardson", atomic_mode="native", _self=cmd)
        assert manager_for(cmd).entries[entry.name] is entry
        reset(cmd)


def test_mesh_budget_rejects_before_allocating_atomic_triangles():
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        sample(cmd)
        entry = apply(cmd)
        with patch(
            "cuemol_style_in_pymol.geometry.sphere",
            side_effect=AssertionError("Unexpected allocation"),
        ):
            with pytest.raises(Exception, match="Atomic meshes exceed cache_mb"):
                cuemol_style("ballstick", atomic_mode="mesh", cache_mb=0.01, _self=cmd)
        assert manager_for(cmd).entries[entry.name] is entry
        reset(cmd)


def test_transparent_native_views_share_mode_and_survive_session_reload():
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        sample(cmd)
        cmd.create("other", "sample", zoom=0)
        cmd.set("transparency_mode", 2)
        apply(cmd, selection="sample", name="view_one", transparency=0.4)
        apply(cmd, selection="other", name="view_two", transparency=0.4)
        assert cmd.get_setting_int("transparency_mode") == 3
        saved = cmd.get_session()
        cmd.reinitialize()
        cmd.set_session(saved)
        assert cmd.get_setting_int("transparency_mode") == 3
        cuemol_style("reset", name="view_one", _self=cmd)
        assert cmd.get_setting_int("transparency_mode") == 3
        cuemol_style("reset", name="view_two", _self=cmd)
        assert cmd.get_setting_int("transparency_mode") == 2
        assert sorted(cmd.get_names("objects")) == ["other", "sample"]
