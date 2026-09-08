"""Geometry invariants independent of a PyMOL installation or a GPU."""

from dataclasses import replace

import numpy as np
import pytest


from cuemol_style_in_pymol import geometry, presets, source, picking
from cuemol_style_in_pymol import mesh as mesh_module


def atom(index, name="CA", ss="S", **kwargs):
    base = source.Atom(
        "peptide",
        index + 1,
        name,
        str(index + 1),
        "ALA",
        "A",
        "",
        "C",
        ss,
        "protein",
        (0.8, 0.25, 0.3),
        1.7,
    )
    return replace(base, **kwargs)


def strand(count=8):
    atoms, coords = [], []
    for i in range(count):
        ca = np.array([i * 3.7, 0.0, 0.0])
        for name, offset in (
            ("CA", [0, 0, 0]),
            ("C", [1, 0, 0]),
            ("O", [1, (-1) ** i, 0]),
        ):
            atoms.append(atom(i, name))
            coords.append(ca + offset)
    return atoms, np.asarray(coords)


@pytest.mark.parametrize("kind", ["rectangle", "ellipse", "fancy"])
def test_winding_normals_and_back_color(kind):
    path = np.array([[0, 0, 0], [2, 0.2, 0.1], [4, 0.4, 0.3]])
    mesh = geometry.sweep(
        path,
        1.2,
        0.2,
        np.tile([0, 1, 0], (3, 1)),
        np.tile([0.8, 0.25, 0.3], (3, 1)),
        [0, 1, 2],
        kind,
        16,
        True,
    )
    triangles = mesh.vertices[mesh.faces]
    face_normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    nonzero = np.linalg.norm(face_normals, axis=1) > 1e-7
    assert np.all(
        np.sum(
            face_normals[nonzero] * mesh.normals[mesh.faces[nonzero]].mean(axis=1),
            axis=1,
        )
        > 0
    )
    assert np.allclose(np.linalg.norm(mesh.normals, axis=1), 1, atol=1e-5)
    assert mesh.colors[:, 1].max() > 0.5
    assert np.isfinite(mesh.edge_data()).all()


@pytest.mark.parametrize("profile_name", ["richardson", "fancy_ribbon", "toon1"])
@pytest.mark.parametrize("hint_sign", [-1, 1])
def test_helix_outer_face_keeps_color_when_frame_sign_changes(profile_name, hint_sign):
    rotation = np.array(
        [[1, 0, 0], [0, np.cos(0.71), -np.sin(0.71)], [0, np.sin(0.71), np.cos(0.71)]]
    )
    atoms, coords = [], []
    for i in range(12):
        theta = np.deg2rad(100 * i)
        ca = np.array([2.3 * np.cos(theta), 2.3 * np.sin(theta), 1.5 * i])
        for name, offset in (
            ("CA", [0, 0, 0]),
            ("C", [0.5, 0, 0]),
            ("O", [0.5, 0, hint_sign]),
        ):
            atoms.append(atom(i, name=name, ss="H"))
            coords.append((ca + offset) @ rotation.T)
    profile = presets.PROFILES[profile_name]
    colors = geometry.atom_colors(atoms, "cuemol")
    mesh = geometry.polymer_mesh(
        atoms, np.asarray(coords), colors, profile, "ribbon", "medium"
    )
    local = mesh.vertices @ rotation
    radial = geometry.unit(local * [1, 1, 0])
    facing = np.sum((mesh.normals @ rotation) * radial, axis=1)
    interior = (mesh.owners // 3 >= 2) & (mesh.owners // 3 <= 9)
    outer, inner = interior & (facing > 0.85), interior & (facing < -0.85)
    assert outer.sum() > 50 and inner.sum() > 50
    np.testing.assert_allclose(
        mesh.colors[outer],
        np.tile(presets.CUEMOL_SECONDARY_COLORS["H"], (outer.sum(), 1)),
        atol=1e-6,
    )
    assert np.all(mesh.colors[inner, 2] > 0.9)


def test_alternating_carbonyls_do_not_twist_sheet_and_arrow_points_forward():
    atoms, coords = strand()
    mesh = geometry.polymer_mesh(
        atoms,
        coords,
        geometry.atom_colors(atoms, "keep"),
        presets.resolve("richardson"),
        "ribbon",
        "medium",
    )
    # The planar input must stay planar through alternating peptide directions.
    assert np.max(np.abs(mesh.vertices[:, 2])) <= 0.201
    shoulder = mesh.vertices[(mesh.vertices[:, 0] > 22) & (mesh.vertices[:, 0] < 24)]
    tip = mesh.vertices[mesh.vertices[:, 0] > 25.5]
    assert np.ptp(shoulder[:, 1]) > np.ptp(tip[:, 1]) * 1.5
    assert mesh.owners.max() == 21


def test_chain_breaks_missing_atoms_and_altlocs_do_not_bridge():
    atoms, coords = strand(8)
    coords[12:] += [50, 0, 0]
    # A high-occupancy B alternate must not replace the main conformer.
    atoms.append(replace(atoms[0], alt="B", occupancy=1.0))
    coords = np.vstack((coords, [1000, 1000, 1000]))
    mesh = geometry.polymer_mesh(
        atoms,
        coords,
        geometry.atom_colors(atoms, "keep"),
        presets.resolve("ribbon"),
        "ribbon",
        "medium",
    )
    triangle = mesh.vertices[mesh.faces]
    assert np.linalg.norm(triangle[:, 1] - triangle[:, 0], axis=1).max() < 6
    assert mesh.vertices.max() < 100
    # Missing a whole residue's CA also separates otherwise adjacent residues.
    atoms[6] = replace(atoms[6], name="CB")
    broken = geometry.polymer_mesh(
        atoms,
        coords,
        geometry.atom_colors(atoms, "keep"),
        presets.resolve("ribbon"),
        "ribbon",
        "medium",
    )
    assert not np.any((broken.vertices[:, 0] > 5) & (broken.vertices[:, 0] < 10))


def test_surface_picker_obeys_native_depth_and_original_owner():
    from types import SimpleNamespace

    mesh = geometry.sphere([0, 0, 0], 0.4, [0.8, 0.25, 0.3], 0, 16)
    drawing = SimpleNamespace(
        matrices=(np.eye(4), np.eye(4), [0, 0, 100, 100]),
        pieces=[SimpleNamespace(mesh=mesh, atoms=(atom(4),))],
    )
    assert picking.hit_at([drawing], 50, 50).index == 5
    assert picking.hit_at([drawing], 50, 50, depth=0.1) is None
    assert picking.hit_at([drawing], 99, 99) is None


def test_degenerate_faces_do_not_create_crease_edges():
    vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    mesh = mesh_module.Mesh(
        vertices, [[0, 0, 1]] * 3, [[1, 0, 0]] * 3, [[0, 1, 2], [0, 0, 1]], [0] * 3
    )
    assert len(mesh.edge_data()) == 3


def test_presets_reject_invalid_controls():
    for arguments in (
        ("bad",),
        ("ribbon", "bad"),
        ("ribbon", "auto", "bad"),
        ("ribbon", "auto", "none", float("nan")),
    ):
        with pytest.raises(ValueError):
            presets.resolve(*arguments)


@pytest.mark.parametrize(
    "representation",
    ["ribbon", "cartoon", "tube", "nucleic", "ballstick", "sticks", "cpk", "surface"],
)
def test_default_cuemol_palette_uses_gui_molecular_painting(representation):
    atoms = [
        atom(0, ss="H"),
        atom(1, ss="S"),
        atom(2, ss="L"),
        atom(3, kind="other", element="S"),
        atom(4, kind="nucleic"),
        atom(5, ss="H", element="N"),
        atom(6, kind="other", ss=""),
    ]
    colors = geometry.atom_colors(atoms, "cuemol", representation)
    nitrogen = (
        [240 / 255, 230 / 255, 140 / 255]
        if representation in ("ribbon", "cartoon", "tube", "nucleic")
        else [0, 0, 1]
    )
    np.testing.assert_allclose(
        colors,
        [
            [240 / 255, 230 / 255, 140 / 255],
            [70 / 255, 130 / 255, 180 / 255],
            [1, 250 / 255, 240 / 255],
            [0, 1, 0],
            [1, 1, 0],
            nitrogen,
            [1, 250 / 255, 240 / 255],
        ],
    )
    np.testing.assert_allclose(
        geometry.atom_colors(atoms[:3], "ss"),
        [[1, 127 / 255, 127 / 255], [127 / 255, 1, 127 / 255], [1, 1, 191 / 255]],
    )
    np.testing.assert_allclose(
        geometry.atom_colors(atoms, "keep"), [a.color for a in atoms]
    )
