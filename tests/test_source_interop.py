"""Real PyMOL checks for managed copies and private molecular sources."""

import numpy as np
import pymol2
import pytest

from cuemol_style_in_pymol import cuemol_style
from cuemol_style_in_pymol.controller import manager_for
from cuemol_style_in_pymol.source import read


@pytest.fixture
def cmd():
    with pymol2.PyMOL() as instance:
        instance.cmd.fragment("ala", "sample")
        instance.cmd.show_as("sticks", "sample")
        yield instance.cmd


def masks(cmd, selection):
    result = []
    cmd.iterate(selection, "result.append(reps)", space={"result": result})
    return result


@pytest.mark.parametrize("hide_names", [0, 1])
def test_private_source_enable_disable_and_reset(cmd, hide_names):
    cmd.set_name("sample", "_private_source")
    cmd.set("hide_underscore_names", hide_names)
    before = masks(cmd, "_private_source")
    cuemol_style("ballstick", "_private_source", quality="low", quiet=1, _self=cmd)
    manager = manager_for(cmd)
    manager.maintenance()
    assert list(manager.active_drawings())
    cmd.disable("_private_source")
    manager.maintenance()
    assert not list(manager.active_drawings())
    cmd.enable("_private_source")
    manager.maintenance()
    assert list(manager.active_drawings())
    cuemol_style("reset", name="all", _self=cmd)
    assert masks(cmd, "_private_source") == before


@pytest.mark.parametrize("reload_session", [False, True])
def test_only_owned_chimerax_copies_are_excluded(cmd, reload_session):
    cmd.create("_cxs_display", "sample")
    cmd.create("_cxs_user_source", "sample")
    cmd._pymol.session.chimerax_style_state = {
        "views": {"example": {"objects": ["_cxs_display"]}}
    }
    if reload_session:
        saved_session = cmd.get_session()
        cmd.reinitialize()
        cmd.set_session(saved_session)
    states, saved = read(cmd, "all")
    assert set(states) == {"sample", "_cxs_user_source"}
    assert {key[0] for key in saved} == set(states)
    with pytest.raises(ValueError, match="no molecular objects"):
        read(cmd, "_cxs_display")


def test_display_copy_deletion_does_not_invalidate_source_view(cmd):
    cmd.create("_cxs_display", "sample")
    cmd._pymol.session.chimerax_style_state = {
        "views": {"example": {"objects": ["_cxs_display"]}}
    }
    original_masks = masks(cmd, "sample")
    copy_masks = masks(cmd, "_cxs_display")
    coordinates = cmd.get_coords("sample").copy()
    entry = cuemol_style("ballstick", quality="low", quiet=1, _self=cmd)
    assert {key[0] for key in entry.saved} == {"sample"}
    assert masks(cmd, "_cxs_display") == copy_masks
    cmd.delete("_cxs_display")
    cmd._pymol.session.chimerax_style_state["views"].clear()
    manager = manager_for(cmd)
    manager.maintenance()
    assert entry.name in manager.entries
    cuemol_style("refresh", _self=cmd)
    assert list(manager.active_drawings())
    cuemol_style("reset", name="all", _self=cmd)
    assert masks(cmd, "sample") == original_masks
    np.testing.assert_array_equal(cmd.get_coords("sample"), coordinates)
