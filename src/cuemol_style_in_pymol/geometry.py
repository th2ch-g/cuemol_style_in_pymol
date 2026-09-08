"""Molecular mesh construction inspired by CueMol's documented sections.

Ribbon and tube axes use natural cubic splines with chord-length knots.
Dimensions follow CueMol renderer styles in angstrom units. Cartoon junctions
and solvent-excluded surfaces remain independent approximations.
"""

import colorsys
from copy import deepcopy
from functools import lru_cache
import xml.etree.ElementTree as ET

import numpy as np

from .mesh import Mesh, merge, unit
from .presets import (
    COIL_COLOR,
    CUEMOL_RADII,
    CUEMOL_ELEMENTS,
    CUEMOL_NUCLEIC_COLOR,
    CUEMOL_OTHER_COLOR,
    CUEMOL_SECONDARY_COLORS,
    QUALITIES,
    SECONDARY_COLORS,
)


@lru_cache(maxsize=32)
def section(kind, detail):
    if kind == "rectangle":
        # Duplicate corners to keep face normals sharp.
        points = np.array(
            [[-1, -1], [1, -1], [1, -1], [1, 1], [1, 1], [-1, 1], [-1, 1], [-1, -1]],
            float,
        )
        normals = np.array(
            [[0, -1], [0, -1], [1, 0], [1, 0], [0, 1], [0, 1], [-1, 0], [-1, 0]], float
        )
        return points, normals
    if kind == "fancy":
        # Two circular rails joined by flat faces, using Fancy1's sharp=0.3.
        theta = 0.3 * np.pi
        angle = np.linspace(
            theta - np.pi / 2, 3 * np.pi / 2 - theta, max(8, detail // 2) + 1
        )
        rail = np.c_[(1.1 + 0.2 * np.sin(angle)) / 1.3, np.cos(angle)]
        normal = np.c_[1.3 * np.sin(angle), 0.2 * np.cos(angle)]
        points = np.vstack((rail, rail[-1], -rail[0], -rail, -rail[-1], rail[0]))[::-1]
        normals = np.vstack((normal, [0, -1], [0, -1], -normal, [0, 1], [0, 1]))[::-1]
        return points, unit(normals)
    t = np.linspace(0, 2 * np.pi, detail, endpoint=False)
    points = np.c_[np.cos(t), np.sin(t)]
    return points, points.copy()


def interpolate(points, samples):
    from scipy.interpolate import CubicSpline

    p = np.asarray(points, float)
    if len(p) == 1:
        return p.copy(), np.array([0.0])
    x = np.arange((len(p) - 1) * samples + 1) / samples
    knots = np.r_[
        0, np.cumsum(np.maximum(np.linalg.norm(np.diff(p, axis=0), axis=1), 1e-8))
    ]
    sample = np.interp(x, np.arange(len(p)), knots)
    return CubicSpline(knots, p, bc_type="natural")(sample), x


@lru_cache(maxsize=64)
def smoothing_operator(count, samples, rho=3.0):
    """Fit the same normalized natural-spline curvature penalty as CueMol."""
    from scipy.interpolate import CubicSpline

    nodes = np.linspace(-1, 1, max(4, count))
    basis = CubicSpline(nodes, np.eye(len(nodes)), bc_type="natural")
    design = basis(np.linspace(-1, 1, count))
    # Two-point Gaussian quadrature integrates products of linear second
    # derivatives exactly on each cubic-spline interval.
    mid = (nodes[1:] + nodes[:-1]) / 2
    half = np.diff(nodes) / 2
    abscissa = (mid[:, None] + half[:, None] * np.array([-1, 1]) / np.sqrt(3)).ravel()
    derivative = basis(abscissa, 2)
    penalty = derivative.T @ (np.repeat(half, 2)[:, None] * derivative)
    scale = np.sum(design**2, axis=0).max()
    weight = 10.0**rho * scale / np.diag(penalty).max()
    normal = design.T @ design + weight * penalty
    normal += np.eye(len(nodes)) * scale * (1 + weight) * 10 * np.finfo(float).eps
    coefficients = np.linalg.solve(normal, design.T)
    return basis(np.linspace(-1, 1, (count - 1) * samples + 1)) @ coefficients


def frames(path, hints, tangent=None):
    tangent = unit(np.gradient(path, axis=0) if tangent is None else tangent)
    side = hints - tangent * np.sum(hints * tangent, axis=1, keepdims=True)
    for i in range(len(side)):
        if np.linalg.norm(side[i]) < 1e-7:
            previous = side[i - 1] if i else np.eye(3)[np.argmin(np.abs(tangent[i]))]
            side[i] = previous - tangent[i] * np.dot(previous, tangent[i])
        if i and np.dot(side[i], side[i - 1]) < 0:
            side[i] *= -1
    side = unit(side)
    return side, unit(np.cross(tangent, side))


def sweep(
    path,
    widths,
    thickness,
    hints,
    colors,
    owners,
    kind,
    detail,
    back=False,
    front=None,
    side_color=False,
    frame=None,
):
    if len(path) < 2:
        return Mesh([], [], [], [], [])
    side, up = frames(path, hints) if frame is None else frame
    shape, sn = section(kind, detail)
    k = len(shape)
    width = np.broadcast_to(widths, (len(path),))
    thick = np.broadcast_to(thickness, (len(path),))
    v = (
        path[:, None]
        + side[:, None] * width[:, None, None] * shape[None, :, 0, None]
        + up[:, None] * thick[:, None, None] * shape[None, :, 1, None]
    )
    n = unit(
        side[:, None] * sn[None, :, 0, None] / np.maximum(width[:, None, None], 1e-6)
        + up[:, None] * sn[None, :, 1, None] / np.maximum(thick[:, None, None], 1e-6)
    )
    col = np.repeat(np.asarray(colors)[:, None, :], k, axis=1)
    if back or side_color:
        # Frame signs can depend on preceding strands and loops. Color helix
        # undersides by the physical inward direction without twisting the mesh.
        polarity = (
            np.where(np.sum(up * front, axis=1) < 0, -1, 1)
            if front is not None
            else np.ones(len(path))
        )
        mask = polarity[:, None] * shape[None, :, 1] < -0.1
        if front is not None:
            mask &= np.linalg.norm(front, axis=1)[:, None] > 1e-7
        if side_color:
            mask = np.broadcast_to(np.abs(sn[:, 0]) > 0.5, col.shape[:2])
        original = col[mask]
        value = original.max(axis=-1, keepdims=True)
        saturation = (value - original.min(axis=-1, keepdims=True)) / np.maximum(
            value, 1e-8
        )
        scale = np.maximum(saturation - 0.4, 0) / np.maximum(saturation, 1e-8)
        col[mask] = value - (value - original) * scale
    j = np.flatnonzero(
        np.linalg.norm(shape - np.roll(shape, -1, axis=0), axis=1) > 1e-8
    )
    a = (np.arange(len(path) - 1)[:, None] * k + j).ravel()
    b = (np.arange(len(path) - 1)[:, None] * k + (j + 1) % k).ravel()
    f = np.concatenate((np.c_[a, b, a + k], np.c_[b, b + k, a + k])).tolist()
    vertices, normals, vertex_colors = (
        v.reshape(-1, 3),
        n.reshape(-1, 3),
        col.reshape(-1, 3),
    )
    ids = np.repeat(owners, k)
    # Independent cap vertices prevent shading the cross-section as a side wall.
    for idx, sign in ((0, -1), (-1, 1)):
        rim = v[idx]
        start = len(vertices)
        normal = unit(np.cross(side[idx], up[idx]))
        vertices = np.vstack((vertices, path[idx], rim))
        normals = np.vstack((normals, np.tile(sign * normal, (k + 1, 1))))
        vertex_colors = np.vstack((vertex_colors, np.tile(colors[idx], (k + 1, 1))))
        ids = np.r_[ids, np.full(k + 1, owners[idx])]
        for j in range(k):
            tri = (start, start + 1 + j, start + 1 + (j + 1) % k)
            f.append(tri if sign > 0 else tri[::-1])
    faces = np.asarray(f)
    # Correct winding against analytic normals, including capped section seams.
    fn = np.cross(
        vertices[faces[:, 1]] - vertices[faces[:, 0]],
        vertices[faces[:, 2]] - vertices[faces[:, 0]],
    )
    reverse = np.sum(fn * normals[faces].mean(axis=1), axis=1) < 0
    faces[reverse] = faces[reverse, ::-1]
    return Mesh(vertices, normals, vertex_colors, faces, ids)


@lru_cache(maxsize=8)
def sphere_template(detail):
    # Rings exclude the poles, whose triangles are added separately.
    phi = np.linspace(0, np.pi, detail // 2 + 1)[1:-1]
    theta = np.linspace(0, 2 * np.pi, detail, endpoint=False)
    v = np.c_[
        np.outer(np.sin(phi), np.cos(theta)).ravel(),
        np.outer(np.sin(phi), np.sin(theta)).ravel(),
        np.repeat(np.cos(phi), detail),
    ]
    f = []
    for i in range(len(phi) - 1):
        for j in range(detail):
            a = i * detail + j
            b = i * detail + (j + 1) % detail
            f.extend(((a, b, a + detail), (b, b + detail, a + detail)))
    top, bottom = len(v), len(v) + 1
    for j in range(detail):
        f.extend(
            (
                (top, (j + 1) % detail, j),
                (
                    bottom,
                    (len(phi) - 1) * detail + j,
                    (len(phi) - 1) * detail + (j + 1) % detail,
                ),
            )
        )
    v = np.vstack((v, [0, 0, 1], [0, 0, -1]))
    f = np.asarray(f)
    reverse = (
        np.sum(
            np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
            * v[f].mean(axis=1),
            axis=1,
        )
        < 0
    )
    f[reverse] = f[reverse, ::-1]
    return v, f


def sphere(center, radius, color, owner, detail):
    v, f = sphere_template(detail)
    return Mesh(
        v * radius + center, v, np.tile(color, (len(v), 1)), f, np.full(len(v), owner)
    )


def bond(a, b, radius, colors, owners, detail):
    if np.linalg.norm(b - a) < 1e-8:
        return Mesh([], [], [], [], [])
    t = unit(b - a)
    h = np.eye(3)[np.argmin(np.abs(t))]
    middle = (a + b) / 2
    # Independent half cylinders keep the element-color boundary sharp.
    return merge(
        [
            sweep(
                np.array([start, end]),
                radius,
                radius,
                np.tile(h, (2, 1)),
                [color, color],
                [owner, owner],
                "ellipse",
                detail,
            )
            for start, end, color, owner in (
                (a, middle, colors[0], owners[0]),
                (middle, b, colors[1], owners[1]),
            )
        ]
    )


def atom_colors(atoms, mode, representation="ribbon"):
    colors = np.array([a.color for a in atoms])
    chains = sorted(set((a.segi, a.chain) for a in atoms))
    elements = {
        "C": (0.45, 0.45, 0.45),
        "N": (0.2, 0.3, 0.9),
        "O": (0.9, 0.15, 0.15),
        "S": (0.95, 0.8, 0.15),
        "P": (1, 0.5, 0.1),
        "H": (0.9, 0.9, 0.9),
    }
    for i, a in enumerate(atoms):
        if mode == "cuemol":
            molecular = (
                CUEMOL_NUCLEIC_COLOR
                if a.kind == "nucleic"
                else CUEMOL_SECONDARY_COLORS.get(a.ss, CUEMOL_OTHER_COLOR)
                if a.kind == "protein"
                else CUEMOL_OTHER_COLOR
            )
            if representation in (
                "ribbon",
                "cartoon",
                "tube",
                "nucleic",
            ) and a.kind in ("protein", "nucleic"):
                colors[i] = molecular
            elif a.element.upper() == "C":
                # DefaultCPKColoring inherits the molecule's paint for carbon.
                colors[i] = molecular
            else:
                colors[i] = CUEMOL_ELEMENTS.get(a.element.upper(), (0.7, 0.7, 0.7))
        elif mode == "chain":
            colors[i] = colorsys.hsv_to_rgb(
                chains.index((a.segi, a.chain)) / max(len(chains), 1), 0.55, 0.9
            )
        elif mode == "rainbow":
            colors[i] = colorsys.hsv_to_rgb(
                (1 - i / max(len(atoms) - 1, 1)) * 0.7, 0.75, 0.95
            )
        elif mode == "ss":
            colors[i] = SECONDARY_COLORS.get(a.ss, COIL_COLOR)
        elif mode == "element":
            colors[i] = elements.get(a.element, (0.7, 0.5, 0.8))
    return colors


def helix_outward(ca, hints, secondary):
    """Find the outside of each helix independently of transported frame signs."""
    outward = np.zeros_like(ca, dtype=float)
    tangent = unit(np.gradient(ca, axis=0))
    start = 0
    while start < len(ca):
        end = start + 1
        while end < len(ca) and secondary[end] == secondary[start]:
            end += 1
        if secondary[start] == "H":
            # Signed carbonyls also orient short or locally straight fragments.
            direction = unit(np.cross(tangent[start:end], hints[start:end]))
            if end - start >= 3:
                local_tangent = unit(np.gradient(ca[start:end], axis=0))
                curvature = -np.gradient(local_tangent, axis=0)
                curvature -= local_tangent * np.sum(
                    curvature * local_tangent, axis=1, keepdims=True
                )
                curved = np.linalg.norm(curvature, axis=1) > 1e-7
                direction[curved] = unit(curvature[curved])
            outward[start:end] = direction
        start = end
    return outward


def polymer_mesh(atoms, coords, colors, profile, representation, quality):
    axial, detail = QUALITIES[quality]
    residues = {}
    for i, a in enumerate(atoms):
        if a.kind in ("protein", "nucleic"):
            residue = residues.setdefault((a.segi, a.chain, a.resi), {})
            old = residue.get(a.name)
            # Prefer blank/A alternate locations, then the largest occupancy.
            priority = (a.alt in ("", "A"), a.occupancy)
            if old is None or priority > (
                atoms[old].alt in ("", "A"),
                atoms[old].occupancy,
            ):
                residue[a.name] = i
    segments = []
    current = []
    last_key = None
    for key, r in residues.items():
        pivot = r.get("CA", r.get("P"))
        if pivot is None:
            if current:
                segments.append(current)
            current = []
            last_key = None
            continue
        if current and (
            key[:2] != last_key[:2]
            or np.linalg.norm(coords[pivot] - coords[current[-1][0]])
            > (4.8 if atoms[pivot].kind == "protein" else 9.0)
        ):
            segments.append(current)
            current = []
        current.append((pivot, r))
        last_key = key
    if current:
        segments.append(current)
    meshes = []
    for segment in segments:
        indices = np.array([p for p, _ in segment])
        ca = coords[indices]
        if len(ca) < 2:
            meshes.append(sphere(ca[0], 0.25, colors[indices[0]], indices[0], detail))
            continue
        seq = np.array([atoms[i].ss for i in indices])
        nucleic = atoms[indices[0]].kind == "nucleic"
        hints = []
        for j, (p, r) in enumerate(segment):
            if seq[j] == "S" and "C" in r and "O" in r:
                hint = coords[r["O"]] - coords[r["C"]]
            elif len(ca) >= 3:
                k = np.clip(j, 1, len(ca) - 2)
                hint = np.cross(ca[k] - ca[k - 1], ca[k + 1] - ca[k])
            else:
                hint = np.array([1.0, 0, 0])
            hints.append(hint)
        axis_points = ca.copy()
        if representation == "ribbon" and not nucleic and len(ca) > 2:
            interior = np.flatnonzero(seq[1:-1] == "S") + 1
            axis_points[interior] = (
                ca[interior] * 0.5 + (ca[interior - 1] + ca[interior + 1]) * 0.25
            )
        path, x = interpolate(axis_points, axial)
        aligned, _ = frames(axis_points, np.asarray(hints))
        normal_path, _ = interpolate(axis_points + aligned, axial)
        hi = normal_path - path
        pick = np.clip(np.floor(x + 0.5).astype(int), 0, len(indices) - 1)
        owner = indices[pick]
        ss = seq[pick]
        front = None
        if profile.back and representation == "ribbon":
            front, _ = interpolate(helix_outward(ca, np.asarray(hints), seq), axial)
            front[ss != "H"] = 0
        fancy = profile.section == "fancy" and representation == "ribbon"
        coil = 0.2 if representation == "cartoon" else 0.25 if fancy else 0.35
        width = np.full(len(path), 1.25 if nucleic else coil)
        thickness = np.full(len(path), 0.5 if nucleic else coil)
        if representation not in ("tube", "nucleic") and not nucleic:
            width[ss == "H"] = 1.3 if fancy else 1.2
            width[ss == "S"] = 1.2 if fancy else 1.4
            thickness[np.isin(ss, ["H", "S"])] = 0.2
            for j, secondary in enumerate(seq):
                if secondary == "S" and (j == len(seq) - 1 or seq[j + 1] != "S"):
                    begin, end = max(0, j - 0.5), min(len(seq) - 1, j + 0.5)
                    region = (x >= begin) & (x <= end)
                    t = (x[region] - begin) / max(end - begin, 1e-8)
                    gamma = (
                        1.0
                        if fancy or representation == "cartoon"
                        else 1.2
                        if profile.section == "ellipse"
                        else 2.2
                    )
                    shoulder = (1.2 if fancy else 1.4) * (1.6 if fancy else 1.8)
                    tip = coil if j < len(seq) - 1 else 0.025
                    width[region] = tip + (shoulder - tip) * (1 - t) ** gamma
        sheet_width = width.copy()
        sheet_thickness = thickness.copy()
        if representation == "cartoon":
            start = 0
            while start < len(seq):
                end = start + 1
                while end < len(seq) and seq[end] == seq[start]:
                    end += 1
                if seq[start] in ("H", "S") and end - start >= 3:
                    part = ca[start:end]
                    rho = (
                        1.0
                        if seq[start] == "S" and profile.section == "ellipse"
                        else 3.0
                    )
                    fitted = smoothing_operator(len(part), axial, rho) @ part
                    projected = fitted[::axial]
                    fitted_hints = (
                        smoothing_operator(len(part), axial, rho) @ aligned[start:end]
                    )
                    mask = (x >= max(0, start - 0.5)) & (
                        x <= min(len(ca) - 1, end - 0.5)
                    )
                    sample_x = np.arange(len(fitted)) / axial + start
                    for axis in range(3):
                        path[mask, axis] = np.interp(x[mask], sample_x, fitted[:, axis])
                        hi[mask, axis] = np.interp(
                            x[mask], sample_x, fitted_hints[:, axis]
                        )
                    before = mask & (x < start)
                    after = mask & (x > end - 1)
                    path[before] += (
                        (x[before] - start)[:, None] * (fitted[1] - fitted[0]) * axial
                    )
                    path[after] += (
                        (x[after] - end + 1)[:, None]
                        * (fitted[-1] - fitted[-2])
                        * axial
                    )
                    if seq[start] == "H":
                        radius = np.linalg.norm(part - projected, axis=1).mean() + 0.2
                        width[mask] = thickness[mask] = max(radius, 0.2)
                start = end
        # Sections change only at secondary-structure boundaries.
        kinds = np.where(np.isin(ss, ["H", "S"]), profile.section, "ellipse").astype(
            object
        )
        if profile.section == "fancy":
            kinds[ss == "S"] = "rectangle"
        if representation in ("tube", "nucleic") or nucleic:
            kinds[:] = "ellipse"
        if representation == "cartoon":
            kinds[ss == "H"] = "ellipse"
        # Use one frame at each shared boundary to keep adjoining sections
        # in the same plane, including direct sheet-to-helix transitions.
        tangent = np.gradient(path, axis=0)
        if representation == "cartoon":
            boundaries = np.r_[0, np.flatnonzero(ss[1:] != ss[:-1]) + 1, len(path) - 1]
            for first, last in zip(boundaries[:-1], boundaries[1:]):
                if ss[first] == "H" and last > first:
                    tangent[first] = path[first + 1] - path[first]
                    tangent[last] = path[last] - path[last - 1]
        side, up = frames(path, hi, tangent)
        begin = 0
        while begin < len(path) - 1:
            end = begin + 1
            while (
                end < len(path) - 1
                and kinds[end] == kinds[begin]
                and ss[end] == ss[begin]
            ):
                end += 1
            sl = slice(begin, end + 1)
            local_width, local_thickness = width[sl], thickness[sl]
            if ss[begin] == "S":
                # A neighboring cylinder must not widen the sheet arrow tip.
                local_width, local_thickness = sheet_width[sl], sheet_thickness[sl]
            if representation == "cartoon" and ss[begin] not in ("H", "S"):
                local_width = local_thickness = 0.2
            meshes.append(
                sweep(
                    path[sl],
                    local_width,
                    local_thickness,
                    hi[sl],
                    colors[owner[sl]],
                    owner[sl],
                    str(kinds[begin]),
                    detail,
                    profile.back and kinds[begin] != "ellipse" and ss[begin] == "H",
                    front[sl] if front is not None else None,
                    side_color=fancy and ss[begin] == "S",
                    frame=(side[sl], up[sl]),
                )
            )
            begin = end
    return merge(meshes)


def nucleic_bases(atoms, coords, colors, detail):
    """DefaultNucl base-pair rods, inferred from hydrogen-bond geometry."""
    from scipy.spatial import cKDTree

    residues = {}
    for i, atom in enumerate(atoms):
        if atom.kind == "nucleic" and atom.alt in ("", "A"):
            residues.setdefault((atom.segi, atom.chain, atom.resi), {})[atom.name] = i
    records = []
    for r in residues.values():
        pivot = r.get("P")
        purine = all(name in r for name in ("C4", "C5", "C8"))
        names = ("C4", "C5", "C8") if purine else ("C2", "C4", "C6")
        end = r.get("N1" if purine else "N3")
        if pivot is None or end is None or not all(name in r for name in names):
            continue
        a, b, c = coords[[r[name] for name in names]]
        normal = unit(np.cross(b - a, c - a))
        records.append((r, pivot, end, normal))
    sites, residue_ids = [], []
    edge_names = {"N1", "N2", "O6", "N6", "N3", "N4", "O2", "O4"}
    for j, (r, _, _, _) in enumerate(records):
        for name, index in r.items():
            if name in edge_names:
                sites.append(index)
                residue_ids.append(j)
    counts = {}
    pairs = {
        frozenset(pair)
        for pair in (
            ("N4", "O6"),
            ("N3", "N1"),
            ("O2", "N2"),
            ("O4", "N6"),
            ("N3", "O6"),
            ("O2", "N1"),
        )
    }
    if sites:
        tree = cKDTree(coords[sites])
        for i, j in sorted(tree.query_pairs(3.7)):
            a, b = sites[i], sites[j]
            ri, rj = residue_ids[i], residue_ids[j]
            if ri == rj or frozenset((atoms[a].name, atoms[b].name)) not in pairs:
                continue
            direction = unit(coords[a] - coords[b])
            if max(
                abs(records[ri][3] @ direction), abs(records[rj][3] @ direction)
            ) > np.sin(np.deg2rad(30)):
                continue
            key = tuple(sorted((ri, rj)))
            counts[key] = counts.get(key, 0) + 1
    partners = {}
    for (i, j), count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        if count >= 2 and i not in partners and j not in partners:
            partners[i], partners[j] = j, i
    meshes = []
    for j, (_, pivot, end, _) in enumerate(records):
        partner = partners.get(j)
        if partner is not None:
            if partner < j:
                continue
            target = records[partner][1]
        else:
            target = end
        meshes.append(
            bond(
                coords[pivot],
                coords[target],
                0.5,
                [
                    colors[pivot],
                    colors[target] if partner is not None else colors[pivot],
                ],
                [pivot, target if partner is not None else pivot],
                detail,
            )
        )
        if partner is None:
            meshes.append(sphere(coords[end], 0.5, colors[pivot], pivot, detail))
    return merge(meshes)


def surface_mesh(model, colors, coords, quality):
    """Use PyMOL's SES algorithm in an isolated, headless instance."""
    import pymol2
    from scipy.spatial import cKDTree

    with pymol2.PyMOL() as p:
        c = p.cmd
        model = deepcopy(model)
        for atom in model.atom:
            atom.vdw = CUEMOL_RADII.get(atom.symbol.upper(), 1.7)
        c.load_model(model, "surface_source")
        c.hide("everything")
        c.show("surface")
        c.set("surface_quality", {"low": 0, "medium": 1, "high": 2}[quality])
        c.set("surface_solvent", 0)
        c.set("solvent_radius", 1.4)
        xml = c.get_collada()
    root = ET.fromstring(xml)
    ns = {"c": "http://www.collada.org/2005/11/COLLADASchema"}
    meshes = []
    for m in root.findall(".//c:mesh", ns):
        sources = {}
        for s in m.findall("c:source", ns):
            arr = s.find("c:float_array", ns)
            accessor = s.find(".//c:accessor", ns)
            if arr is not None:
                sources[s.attrib["id"]] = np.fromstring(
                    arr.text or "", sep=" "
                ).reshape(-1, int(accessor.attrib.get("stride", 3)))
        vert = {
            v.attrib["id"]: v.find("c:input", ns).attrib["source"][1:]
            for v in m.findall("c:vertices", ns)
        }
        for prim in list(m):
            if prim.tag.rsplit("}", 1)[-1] not in ("triangles", "polylist"):
                continue
            inp = {
                i.attrib["semantic"]: (
                    i.attrib["source"][1:],
                    int(i.attrib.get("offset", 0)),
                )
                for i in prim.findall("c:input", ns)
            }
            stride = max(off for _, off in inp.values()) + 1
            indices = np.fromstring(
                prim.find("c:p", ns).text, sep=" ", dtype=int
            ).reshape(-1, stride)
            source, offset = inp["VERTEX"]
            v = sources[vert.get(source, source)][indices[:, offset], :3]
            counts = prim.find("c:vcount", ns)
            counts = (
                np.fromstring(counts.text, sep=" ", dtype=int)
                if counts is not None
                else np.full(len(v) // 3, 3)
            )
            f = []
            first = 0
            for count in counts:
                f.extend((first, first + j, first + j + 1) for j in range(1, count - 1))
                first += count
            f = np.asarray(f)
            if "NORMAL" in inp:
                source, offset = inp["NORMAL"]
                normals = sources[source][indices[:, offset], :3]
            else:
                normals = np.zeros_like(v)
                np.add.at(
                    normals,
                    f.ravel(),
                    np.repeat(
                        np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]),
                        3,
                        axis=0,
                    ),
                )
                normals = unit(normals)
            owners = cKDTree(coords).query(v)[1]
            meshes.append(Mesh(v, normals, colors[owners], f, owners))
    return merge(meshes)


def build(atoms, coords, bonds, model, profile, representation, quality, color_mode):
    colors = atom_colors(atoms, color_mode, representation)
    detail = QUALITIES[quality][1]
    if representation == "surface":
        return surface_mesh(model, colors, coords, quality)
    parts = []
    if representation in ("ribbon", "cartoon", "tube", "nucleic"):
        parts.append(
            polymer_mesh(atoms, coords, colors, profile, representation, quality)
        )
        parts.append(nucleic_bases(atoms, coords, colors, detail))
        chosen = {i for i, a in enumerate(atoms) if a.kind == "other"}
    else:
        chosen = set(range(len(atoms)))
    for i in sorted(chosen):
        a = atoms[i]
        if representation == "cpk":
            radius = CUEMOL_RADII.get(a.element.upper(), 1.7)
        elif representation == "sticks":
            radius = 0.2
        else:
            radius = 0.3
        parts.append(sphere(coords[i], radius, colors[i], i, detail))
    if representation != "cpk":
        for i, j in bonds:
            if i in chosen and j in chosen:
                parts.append(
                    bond(
                        coords[i],
                        coords[j],
                        0.2,
                        [colors[i], colors[j]],
                        [i, j],
                        detail,
                    )
                )
    return merge(parts)
