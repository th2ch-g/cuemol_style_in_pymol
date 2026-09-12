"""CueMol sections, chord-length ribbon splines, and standalone EDTSurf meshes."""

import colorsys
import itertools
from functools import lru_cache

import numpy as np

from .mesh import Mesh, merge, unit
from .presets import (
    COIL_COLOR,
    CUEMOL_ELEMENTS,
    CUEMOL_NUCLEIC_COLOR,
    CUEMOL_OTHER_COLOR,
    CUEMOL_RADII,
    CUEMOL_SECONDARY_COLORS,
    QUALITIES,
    SECONDARY_COLORS,
)


@lru_cache(maxsize=64)
def section(kind, detail, ratio=6.5):
    """CueMol cross sections in normalized (width, thickness) coordinates."""
    if kind == "rectangle":
        count = max(1, int(detail * ratio / (2 * (ratio + 1))))
        points, normals = [], []
        for start, end, normal in (
            ((-1, -1), (1, -1), (0, -1)),
            ((1, -1), (1, 1), (1, 0)),
            ((1, 1), (-1, 1), (0, 1)),
            ((-1, 1), (-1, -1), (-1, 0)),
        ):
            size = count if start[1] == end[1] else max(1, int(count / ratio))
            points.extend(np.linspace(start, end, size + 1))
            normals.extend([normal] * (size + 1))
        return np.asarray(points), np.asarray(normals)
    if kind == "fancy":
        theta = 0.3 * np.pi
        by = max(ratio - 1 - np.cos(theta), 0)
        cy = by + np.cos(theta)
        arc_length = 2 * (np.pi - theta)
        spacing = (4 * by + 2 * arc_length) / detail
        arc_steps = max(4, int(arc_length / spacing)) * 2
        flat_steps = max(1, int(2 * by / spacing))
        angle = np.linspace(theta - np.pi / 2, 3 * np.pi / 2 - theta, arc_steps + 1)
        rail = np.c_[(np.sin(angle) + cy) / ratio, np.cos(angle)]
        normal = np.c_[ratio * np.sin(angle), np.cos(angle)]
        front = np.c_[
            np.linspace(by, -by, flat_steps + 1) / ratio,
            np.full(flat_steps + 1, -np.sin(theta)),
        ]
        points = np.vstack((rail, front, -rail, -front))
        normals = np.vstack(
            (
                normal,
                np.tile([0, -1], (len(front), 1)),
                -normal,
                np.tile([0, 1], (len(front), 1)),
            )
        )
        return points, unit(normals)
    phi = np.linspace(0, 2 * np.pi, detail + 1)
    angles = np.unwrap(np.arctan2(ratio * np.sin(phi), np.cos(phi)))
    samples = []
    for start, end in itertools.pairwise(angles):
        divisions = max(1, int(np.ceil((end - start) / (2 * np.pi / detail + 0.001))))
        samples.extend(np.linspace(start, end, divisions, endpoint=False))
    points = np.c_[np.sin(samples), np.cos(samples)]
    return points, points.copy()


def interpolate(points, samples, derivative=0):
    from scipy.interpolate import CubicSpline

    p = np.asarray(points, float)
    if len(p) == 1:
        return p.copy(), np.array([0.0])
    x = np.arange((len(p) - 1) * samples + 1) / samples
    knots = np.r_[
        0, np.cumsum(np.maximum(np.linalg.norm(np.diff(p, axis=0), axis=1), 1e-8))
    ]
    sample = np.interp(x, np.arange(len(p)), knots)
    result = CubicSpline(knots, p, bc_type="natural")(sample, derivative)
    if derivative:
        interval = np.minimum(np.floor(x).astype(int), len(p) - 2)
        result *= np.diff(knots)[interval, None] ** derivative
    return result, x


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


def transport_hints(path, hints, tangent=None):
    tangent = unit(np.gradient(path, axis=0) if tangent is None else tangent)
    aligned = unit(np.asarray(hints, dtype=float)).copy()
    for i in range(len(aligned)):
        if np.linalg.norm(aligned[i]) < 1e-7:
            aligned[i] = (
                aligned[i - 1] if i else np.eye(3)[np.argmin(np.abs(tangent[i]))]
            )
        if i:
            axis = np.cross(tangent[i - 1], tangent[i])
            cosine = np.clip(tangent[i - 1] @ tangent[i], -1, 1)
            previous = aligned[i - 1]
            transported = previous + np.cross(axis, previous)
            if cosine > -0.999999:
                transported += np.cross(axis, np.cross(axis, previous)) / (1 + cosine)
            if np.dot(aligned[i], transported) < 0:
                aligned[i] *= -1
    return aligned


def frames(path, hints, tangent=None, ribbon=False):
    tangent = unit(np.gradient(path, axis=0) if tangent is None else tangent)
    if ribbon:
        side = np.asarray(hints).copy()
    else:
        side = hints - tangent * np.sum(hints * tangent, axis=1, keepdims=True)
        side = transport_hints(path, side, tangent)
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
    caps=(True, True),
):
    if len(path) < 2:
        return Mesh([], [], [], [], [])
    side, up = frames(path, hints) if frame is None else frame
    width = np.broadcast_to(widths, (len(path),))
    thick = np.broadcast_to(thickness, (len(path),))
    ratio = (
        6.5 if kind == "fancy" else float(np.median(width / np.maximum(thick, 1e-6)))
    )
    shape, sn = section(kind, detail, round(ratio, 6))
    k = len(shape)
    v = (
        path[:, None]
        + side[:, None] * width[:, None, None] * shape[None, :, 0, None]
        + up[:, None] * thick[:, None, None] * shape[None, :, 1, None]
    )
    n = unit(
        side[:, None] * sn[None, :, 0, None] / np.maximum(width[:, None, None], 1e-6)
        + up[:, None] * sn[None, :, 1, None] / np.maximum(thick[:, None, None], 1e-6)
    )
    velocity = np.gradient(path, axis=0)
    speed2 = np.maximum(np.sum(velocity**2, axis=1), 1e-12)
    physical_sn = unit(
        sn / [max(float(np.median(width)), 1e-6), max(float(np.median(thick)), 1e-6)]
    )
    slope = -(
        np.gradient(width)[:, None] * np.abs(physical_sn[None, :, 0])
        + np.gradient(thick)[:, None] * np.abs(physical_sn[None, :, 1])
    )
    n = unit(n + velocity[:, None] * (slope / speed2[:, None])[:, :, None])
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
        if kind == "fancy":
            mask &= np.abs(sn[None, :, 1]) > 0.999
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
    for enabled, (idx, sign) in zip(caps, ((0, -1), (-1, 1))):
        if not enabled:
            continue
        rim = v[idx]
        start = len(vertices)
        normal = unit(np.cross(side[idx], up[idx]))
        if enabled == "sphere":
            extent = max(np.linalg.norm(rim[0] - path[idx]), 1e-8)
            t = np.linspace(0, 1, 6)
            radius = np.sqrt(1 - t * t)
            cap = (
                path[idx]
                + radius[:, None, None] * (rim - path[idx])
                + sign * extent * t[:, None, None] * normal
            )
            cap_normals = unit(
                radius[:, None, None] * n[idx] + sign * t[:, None, None] * normal
            )
            vertices = np.vstack((vertices, cap.reshape(-1, 3)))
            normals = np.vstack((normals, cap_normals.reshape(-1, 3)))
            vertex_colors = np.vstack((vertex_colors, np.tile(colors[idx], (6 * k, 1))))
            ids = np.r_[ids, np.full(6 * k, owners[idx])]
            for ring in range(5):
                for j in range(k):
                    a, b = start + ring * k + j, start + ring * k + (j + 1) % k
                    f.extend(((a, b, a + k), (b, b + k, a + k)))
            continue
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
    v, f = sphere_template(detail * 3)
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
                detail * 3,
            )
            for start, end, color, owner in (
                (a, middle, colors[0], owners[0]),
                (middle, b, colors[1], owners[1]),
            )
        ]
    )


def atom_colors(atoms, mode, representation="ribbon"):
    colors = np.array([a.color for a in atoms])
    chains = sorted({(a.segi, a.chain) for a in atoms})
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
        if representation == "cartoon" and atoms[indices[0]].kind == "protein":
            from .cartoon import build

            meshes.extend(build(segment, atoms, coords, colors, profile, axial, detail))
            continue
        if (
            representation == "ribbon"
            and atoms[indices[0]].kind == "protein"
            and len(ca) >= 3
        ):
            from .ribbon import build

            meshes.extend(build(segment, atoms, coords, colors, profile, axial, detail))
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
        axis_tangent, _ = interpolate(axis_points, axial, 1)
        aligned = transport_hints(axis_points, np.asarray(hints), axis_tangent[::axial])
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
                first = j == 0 or seq[j - 1] != secondary
                last = j == len(seq) - 1 or seq[j + 1] != secondary
                if representation == "ribbon" and secondary in ("H", "S"):
                    begin, end = max(0, j - 0.5), min(len(seq) - 1, j + 0.5)
                    region = (x >= begin) & (x <= end)
                    t = (x[region] - begin) / max(end - begin, 1e-8)
                    blend = np.where(
                        t <= 0.5, (2 * t) ** 2.2 / 2, 1 - (2 * (1 - t)) ** 2.2 / 2
                    )
                    if first or (last and secondary == "H"):
                        weight = blend if first else 1 - blend
                        full_width = (
                            (1.3 if fancy else 1.2)
                            if secondary == "H"
                            else (1.2 if fancy else 1.4)
                        )
                        width[region] = coil + (full_width - coil) * weight
                        thickness[region] = coil + (0.2 - coil) * weight
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
                    tip = coil
                    width[region] = tip + (shoulder - tip) * (1 - t) ** gamma
                    thickness[region] = 0.2 + (coil - 0.2) * np.where(
                        t <= 0.5, (2 * t) ** gamma / 2, 1 - (2 * (1 - t)) ** gamma / 2
                    )
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
        tangent = (
            axis_tangent if representation != "cartoon" else np.gradient(path, axis=0)
        )
        if representation == "cartoon":
            boundaries = np.r_[0, np.flatnonzero(ss[1:] != ss[:-1]) + 1, len(path) - 1]
            for first, last in itertools.pairwise(boundaries):
                if ss[first] == "H" and last > first:
                    tangent[first] = path[first + 1] - path[first]
                    tangent[last] = path[last] - path[last - 1]
        side, up = frames(
            path, hi, tangent, ribbon=representation == "ribbon" and not nucleic
        )
        sampled_colors = colors[owner].copy()
        if representation != "cartoon":
            for channel in range(3):
                sampled_colors[:, channel] = np.interp(
                    x, np.arange(len(indices)), colors[indices, channel]
                )
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
                    sampled_colors[sl],
                    owner[sl],
                    str(kinds[begin]),
                    detail,
                    profile.back and kinds[begin] != "ellipse" and ss[begin] == "H",
                    front[sl] if front is not None else None,
                    side_color=fancy and ss[begin] == "S",
                    frame=(side[sl], up[sl]),
                    caps=(
                        "sphere" if begin == 0 else False,
                        "sphere" if end == len(path) - 1 else False,
                    ),
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
        for i, neighbors in enumerate(tree.query_ball_point(coords[sites], 3.7)):
            a, ri = sites[i], residue_ids[i]
            candidates = []
            for j in neighbors:
                b, rj = sites[j], residue_ids[j]
                if ri == rj or frozenset((atoms[a].name, atoms[b].name)) not in pairs:
                    continue
                delta = coords[a] - coords[b]
                direction = unit(delta)
                if max(
                    abs(records[ri][3] @ direction), abs(records[rj][3] @ direction)
                ) <= np.sin(np.deg2rad(30)):
                    candidates.append((np.dot(delta, delta), rj))
            if candidates:
                _, partner = min(candidates)
                counts[ri, partner] = counts.get((ri, partner), 0) + 1
    partners = {}
    for i in range(len(records)):
        candidates = [(count, -j) for (row, j), count in counts.items() if row == i]
        if not candidates:
            continue
        count, partner = max(candidates)
        j = -partner
        if 2 <= count <= 3 and i not in partners and j not in partners:
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
    """Generate CueMol's EDTSurf solvent-excluded surface and voxel owners."""
    from ._edtsurf import surface

    elements = {name: i for i, name in enumerate(("H", "C", "N", "O", "S", "P"))}
    types = [elements.get(atom.symbol.upper(), 6) for atom in model.atom]
    vertices, normals, faces, owners = surface(
        coords, types, {"low": 3, "medium": 6, "high": 10}[quality]
    )
    return Mesh(vertices, normals, colors[owners], faces, owners)


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
