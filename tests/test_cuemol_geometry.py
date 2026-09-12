"""Geometry invariants independent of a PyMOL installation or a GPU."""

from dataclasses import replace

import numpy as np
import pytest

from cuemol_style_in_pymol import geometry, picking, presets, source
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
    if profile.back:
        # Fancy rails retain their pigment; only the flat underside is lightened.
        assert np.count_nonzero(mesh.colors[inner, 2] > 0.9) > inner.sum() / 2
    else:
        np.testing.assert_allclose(
            mesh.colors[inner], np.tile(colors[0], (inner.sum(), 1))
        )


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
    interior = (mesh.vertices[:, 0] > 4) & (mesh.vertices[:, 0] < 22)
    assert np.max(np.abs(mesh.vertices[interior, 2])) <= 0.201
    assert np.max(np.abs(mesh.vertices[:, 2])) <= 0.251
    shoulder = mesh.vertices[(mesh.vertices[:, 0] > 24) & (mesh.vertices[:, 0] < 24.1)]
    tip = mesh.vertices[mesh.vertices[:, 0] > 27.7]
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


def test_natural_chord_spline_reference_midpoint():
    points, _ = geometry.interpolate([[0, 0, 0], [1, 1, 0], [2, 0, 0]], 2)
    np.testing.assert_allclose(points[1], [0.5, 0.6875, 0], atol=1e-12)
    np.testing.assert_allclose(points[3], [1.5, 0.6875, 0], atol=1e-12)


def test_penalized_axis_preserves_lines_and_smooths_helix():
    t = np.linspace(0, 1, 12)
    line = np.c_[t, 2 * t, -3 * t]
    op = geometry.smoothing_operator(12, 8)
    fit = op @ line
    np.testing.assert_allclose(fit[::8], line, atol=1e-8)
    wave = line + np.c_[np.sin(t * 6 * np.pi), np.zeros((12, 2))]
    assert np.linalg.norm((op @ wave)[::8] - line) < np.linalg.norm(wave - line) / 2


@pytest.mark.parametrize(
    "style,width", [("ribbon", 1.4), ("round_ribbon", 1.4), ("fancy_ribbon", 1.2)]
)
def test_sheet_body_width_matches_cuemol(style, width):
    atoms, coords = strand()
    mesh = geometry.polymer_mesh(
        atoms,
        coords,
        geometry.atom_colors(atoms, "keep"),
        presets.PROFILES[style],
        "ribbon",
        "medium",
    )
    body = (mesh.vertices[:, 0] > 5) & (mesh.vertices[:, 0] < 15)
    np.testing.assert_allclose(np.abs(mesh.vertices[body, 1]).max(), width, atol=1e-6)
    if style == "fancy_ribbon":
        # Fancy sheets have flat faces and desaturated side walls, not oval sections.
        normals = mesh.normals[body]
        assert np.all(
            np.isclose(np.abs(normals[:, 1]), 1) | np.isclose(np.abs(normals[:, 2]), 1)
        )
        front = body & (mesh.normals[:, 2] > 0.99)
        np.testing.assert_allclose(
            mesh.colors[front], np.tile(atoms[0].color, (front.sum(), 1)), atol=1e-6
        )
        assert mesh.colors[body & (mesh.normals[:, 1] > 0.99), 1].mean() > 0.4


def test_fancy_section_has_circular_rails_and_flat_middle():
    points, _ = geometry.section("fancy", 24)
    physical = points * [1.3, 0.2]
    np.testing.assert_allclose(np.abs(physical[:, 0]).max(), 1.3, atol=1e-12)
    inner = np.abs(physical[:, 0]) < 1.0
    np.testing.assert_allclose(
        np.abs(physical[inner, 1]), 0.2 * np.sin(0.3 * np.pi), atol=1e-12
    )


def test_tube_and_nucleic_backbone_dimensions():
    atoms = [atom(i, ss="") for i in range(4)]
    coords = np.c_[np.arange(4) * 3.7, np.zeros((4, 2))]
    colors = geometry.atom_colors(atoms, "keep")
    tube = geometry.polymer_mesh(
        atoms, coords, colors, presets.PROFILES["tube"], "tube", "medium"
    )
    np.testing.assert_allclose(
        np.abs(tube.vertices[:, 1:]).max(axis=0), [0.35, 0.35], atol=1e-6
    )
    nucleic = [replace(a, name="P", kind="nucleic") for a in atoms]
    backbone = geometry.polymer_mesh(
        nucleic, coords, colors, presets.PROFILES["nucleic"], "nucleic", "medium"
    )
    np.testing.assert_allclose(
        np.sort(np.abs(backbone.vertices[:, 1:]).max(axis=0)), [0.5, 1.25], atol=1e-6
    )


def test_atomic_radii_and_sharp_bond_colors():
    atoms = [atom(0, kind="other", element="H", vdw=9.0)]
    coords = np.zeros((1, 3))
    for representation, radius in (("ballstick", 0.3), ("sticks", 0.2), ("cpk", 1.2)):
        mesh = geometry.build(
            atoms,
            coords,
            [],
            None,
            presets.PROFILES["default"],
            representation,
            "medium",
            "keep",
        )
        np.testing.assert_allclose(
            np.linalg.norm(mesh.vertices, axis=1), radius, atol=1e-6
        )
    bond = geometry.bond(
        np.array([0, 0, 0]),
        np.array([4, 0, 0]),
        0.2,
        [[1, 0, 0], [0, 0, 1]],
        [0, 1],
        16,
    )
    np.testing.assert_allclose(
        bond.colors[bond.vertices[:, 0] < 2],
        np.tile([1, 0, 0], ((bond.vertices[:, 0] < 2).sum(), 1)),
    )
    np.testing.assert_allclose(
        bond.colors[bond.vertices[:, 0] > 2],
        np.tile([0, 0, 1], ((bond.vertices[:, 0] > 2).sum(), 1)),
    )


def test_richardson_native_average_retains_paper_and_pigment():
    from cuemol_style_in_pymol.materials import drawing_tone, pencil_average

    tone = drawing_tone(np.array([[0, 0, 1], [1, 0, 0], [-1, 0, 0]]))
    assert tone[0] > 0.92
    assert tone[2] < tone[1] < tone[0]
    colors = pencil_average(np.tile([0.2, 0.4, 0.8], (3, 1)), tone)
    np.testing.assert_allclose(colors[0], np.array([240, 236, 221]) / 255)
    assert colors[2].mean() < colors[0].mean()
    assert colors[2, 2] > colors[2, 0]


@pytest.mark.parametrize(
    "style,representation",
    [
        ("cartoon", "cartoon"),
        ("round_cartoon", "cartoon"),
    ],
)
def test_direct_sheet_helix_junction_has_a_shared_plane_and_narrow_tip(
    monkeypatch, style, representation
):
    atoms, coords = [], []
    for i in range(15):
        theta = np.deg2rad(100 * (i - 5))
        ca = (
            np.array([2.3 - 3.7 * (5 - i), 0.0, 0.0])
            if i < 5
            else np.array([2.3 * np.cos(theta), 2.3 * np.sin(theta), 1.5 * (i - 5)])
        )
        for name, offset in (("CA", [0, 0, 0]), ("C", [1, 0, 0]), ("O", [1, 1, 0])):
            atoms.append(atom(i, name=name, ss="S" if i < 5 else "H"))
            coords.append(ca + offset)
    sections = []
    original = geometry.sweep

    def capture(*args, **kwargs):
        mesh = original(*args, **kwargs)
        ratio = (
            6.5
            if args[6] == "fancy"
            else float(np.median(np.asarray(args[1]) / np.asarray(args[2])))
        )
        count = len(geometry.section(args[6], args[7], round(ratio, 6))[0])
        sections.append((mesh, count, np.asarray(args[0])))
        return mesh

    monkeypatch.setattr(geometry, "sweep", capture)
    geometry.polymer_mesh(
        atoms,
        np.asarray(coords),
        geometry.atom_colors(atoms, "cuemol"),
        presets.PROFILES[style],
        representation,
        "medium",
    )
    assert len(sections) == 2
    sheet, sk, sheet_axis = sections[0]
    helix, hk, helix_axis = sections[1]
    sheet_center = sheet_axis[-1]
    sheet_rim = sheet.vertices[(len(sheet_axis) - 1) * sk : len(sheet_axis) * sk]
    helix_center = helix_axis[0]
    helix_normal = geometry.unit(
        np.cross(
            helix.vertices[1] - helix.vertices[0],
            helix.vertices[hk // 2] - helix.vertices[0],
        )
    )
    if representation == "cartoon":
        # Ribbon2Renderer fits each element independently, including flanks.
        radius = np.linalg.norm(helix.vertices[:hk] - helix_center, axis=1).min()
        assert np.linalg.norm(sheet_rim - helix_center, axis=1).max() < radius
    else:
        np.testing.assert_allclose(sheet_center, helix_center, atol=1e-6)
        assert np.max(np.abs((sheet_rim - helix_center) @ helix_normal)) < 1e-5
    assert np.linalg.norm(sheet_rim - sheet_center, axis=1).max() < 0.5


@pytest.mark.parametrize("style", ["ribbon", "round_ribbon", "fancy_ribbon"])
def test_ribbon_arrow_has_a_flat_shoulder_and_narrow_tip(style):
    atoms, coords = strand(8)
    mesh = geometry.polymer_mesh(
        atoms,
        coords,
        geometry.atom_colors(atoms, "keep"),
        presets.resolve(style),
        "ribbon",
        "medium",
    )
    triangles = mesh.vertices[mesh.faces]
    shoulder = np.all(np.isclose(triangles[:, :, 0], 6.5 * 3.7, atol=1e-5), axis=1)
    area = np.linalg.norm(
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]),
        axis=1,
    )
    assert area[shoulder].sum() > 0.05
    width = 1.2 * 1.6 if style == "fancy_ribbon" else 1.4 * 1.8
    assert np.max(np.abs(triangles[shoulder, :, 1])) == pytest.approx(width, abs=1e-5)
    boundary = mesh.vertices[np.isclose(mesh.vertices[:, 0], 7.5 * 3.7, atol=1e-5)]
    assert len(boundary) > 10
    coil = 0.25 if style == "fancy_ribbon" else 0.35
    assert np.max(np.abs(boundary[:, 1:])) <= coil + 1e-5
