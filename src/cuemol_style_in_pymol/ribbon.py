"""Residue-local ribbon junctions and section normals."""

import numpy as np
from scipy.interpolate import CubicSpline

from .mesh import Mesh, unit


def chord_spline(points):
    knots = np.r_[
        0, np.cumsum(np.maximum(np.linalg.norm(np.diff(points, axis=0), axis=1), 1e-8))
    ]
    spline = CubicSpline(knots, points, bc_type="natural")

    def evaluate(parameters, derivative=0):
        indices = np.clip(np.floor(parameters).astype(int), 0, len(points) - 2)
        lengths = np.diff(knots)[indices]
        t = knots[indices] + (parameters - indices) * lengths
        return spline(t, derivative) * lengths[:, None] ** derivative

    return evaluate


def junction(axial, width, coil, *, head=False, arrow=False, gamma=2.2):
    """Return parameter, physical half-axes, and their analytic derivatives."""
    t = np.linspace(0, 1, axial + 1)
    rho = np.where(t <= 0.5, (2 * t) ** gamma / 2, 1 - (2 * (1 - t)) ** gamma / 2)
    drho = gamma * (2 * np.minimum(t, 1 - t)) ** (gamma - 1)
    if arrow:
        shoulder = width * (1.6 if gamma == 1 else 1.8)
        w = coil + (shoulder - coil) * (1 - t) ** gamma
        dw = -(shoulder - coil) * gamma * np.maximum(1 - t, 1e-4) ** (gamma - 1)
        h = 0.2 + (coil - 0.2) * rho
        dh = (coil - 0.2) * drho
        return (
            np.r_[0, t],
            np.c_[np.r_[0.2, h], np.r_[width, w]],
            np.c_[np.r_[0, dh], np.r_[0, dw]],
        )
    start, finish = (
        ((0.2, width), (coil, coil)) if head else ((coil, coil), (0.2, width))
    )
    delta = np.subtract(finish, start)
    return t, np.asarray(start) + rho[:, None] * delta, drho[:, None] * delta


def section(kind, width, thickness, detail):
    from .geometry import section as normalized_section

    ratio = width / thickness
    if kind == "rectangle":
        nx = max(1, int(detail / (2 * (ratio + 1))))
        ny = max(1, int(detail * ratio / (2 * (ratio + 1))))
        points, normals, back, side = [], [], [], []
        for a, b, n, count, is_back, is_side in (
            ((1, 1), (-1, 1), (0, 1), nx, False, True),
            ((-1, 1), (-1, -1), (-1, 0), ny, False, False),
            ((-1, -1), (1, -1), (0, -1), nx, False, True),
            ((1, -1), (1, 1), (1, 0), ny, True, False),
        ):
            points.extend(np.linspace(a, b, count + 1) * [thickness, width])
            normals.extend([n] * (count + 1))
            back.extend([is_back] * (count + 1))
            side.extend([is_side] * (count + 1))
        return tuple(np.asarray(a) for a in (points, normals, back, side))
    shape, _ = normalized_section(kind, detail, ratio)
    points = shape[:, ::-1] * [thickness, width]
    if kind == "fancy":
        theta = 0.3 * np.pi
        cy = ratio - 1
        normals = unit(
            np.c_[
                points[:, 0] / thickness,
                points[:, 1] / thickness - np.sign(points[:, 1]) * cy,
            ]
        )
        # The section table duplicates the rail/flat corners with distinct normals.
        spacing = (4 * (cy - np.cos(theta)) + 4 * (np.pi - theta)) / detail
        arc = max(4, int(2 * (np.pi - theta) / spacing)) * 2 + 1
        line = max(1, int(2 * (cy - np.cos(theta)) / spacing)) + 1
        flat = np.zeros(len(points), bool)
        flat[arc : arc + line] = True
        flat[2 * arc + line :] = True
        normals[flat] = np.c_[np.sign(points[flat, 0]), np.zeros(flat.sum())]
        back = flat & (points[:, 0] > 0)
        return points, normals, back, np.zeros(len(points), bool)
    normals = unit(points / [thickness**2, width**2])
    return points, normals, np.zeros(len(points), bool), np.zeros(len(points), bool)


def build(segment, atoms, coords, colors, profile, axial, detail):
    from .geometry import transport_hints

    indices = np.array([p for p, _ in segment])
    ca = coords[indices]
    secondary = np.array([atoms[i].ss for i in indices])
    points = ca.copy()
    interior = np.flatnonzero(secondary[1:-1] == "S") + 1
    points[interior] = ca[interior] * 0.5 + (ca[interior - 1] + ca[interior + 1]) * 0.25
    axis = chord_spline(points)
    hints = []
    for j, (_, residue) in enumerate(segment):
        if secondary[j] == "S" and "C" in residue and "O" in residue:
            hint = coords[residue["O"]] - coords[residue["C"]]
        elif len(ca) > 2:
            k = np.clip(j, 1, len(ca) - 2)
            hint = np.cross(ca[k] - ca[k - 1], ca[k + 1] - ca[k])
        else:
            hint = np.array([1.0, 0.0, 0.0])
        hints.append(hint if np.linalg.norm(hint) >= 1e-4 else [1.0, 0.0, 0.0])
    hints = unit(hints)
    aligned = transport_hints(points, hints, axis(np.arange(len(points)), 1))
    inverted = np.sum(hints * aligned, axis=1) < 0
    normal = chord_spline(points + aligned)
    fancy = profile.section == "fancy"
    coil = 0.25 if fancy else 0.35
    meshes = []

    for j, secondary_type in enumerate(secondary):
        first = j == 0 or secondary[j - 1] != secondary_type
        last = j == len(ca) - 1 or secondary[j + 1] != secondary_type
        structured = secondary_type in ("H", "S")
        kind = (
            ("rectangle" if fancy and secondary_type == "S" else profile.section)
            if structured
            else "ellipse"
        )
        width = (
            (
                (1.3 if fancy else 1.2)
                if secondary_type == "H"
                else (1.2 if fancy else 1.4)
            )
            if structured
            else coil
        )
        thickness = 0.2 if structured else coil
        start, end = max(0, j - 0.5), min(len(ca) - 1, j + 0.5)
        arrow = secondary_type == "S" and last and not first
        if structured and (first or last):
            gamma = (
                (1.0 if fancy else 1.2 if kind == "ellipse" else 2.2) if arrow else 2.2
            )
            t, sizes, slopes = junction(
                axial, width, coil, head=not first, arrow=arrow, gamma=gamma
            )
            x = start + t
        else:
            x = np.linspace(start, end, max(1, int((end - start) * axial)) + 1)
            sizes = np.tile([thickness, width], (len(x), 1))
            slopes = np.zeros_like(sizes)
        path, tangent = axis(x), axis(x, 1)
        e2 = normal(x) - path
        e1 = unit(np.cross(e2, tangent))
        cross, normals, backs, sides = section(kind, width, thickness, detail)
        scale, derivative = sizes / [thickness, width], slopes / [thickness, width]
        vertices = (
            path[:, None]
            + e1[:, None] * (cross[None, :, 0] * scale[:, None, 0])[:, :, None]
            + e2[:, None] * (cross[None, :, 1] * scale[:, None, 1])[:, :, None]
        )
        n = (
            e1[:, None] * (normals[None, :, 0] * scale[:, None, 1])[:, :, None]
            + e2[:, None] * (normals[None, :, 1] * scale[:, None, 0])[:, :, None]
        )
        dn = (
            normals[None, :, 0]
            * scale[:, None, 1]
            * cross[None, :, 0]
            * derivative[:, None, 0]
            - normals[None, :, 1]
            * scale[:, None, 0]
            * cross[None, :, 1]
            * derivative[:, None, 1]
        )
        n = unit(
            n
            + tangent[:, None]
            * (dn / np.maximum(np.sum(tangent**2, axis=1), 1e-12)[:, None])[:, :, None]
        )
        pigment = np.array(
            [
                np.interp(x, np.arange(len(ca)), colors[indices, channel])
                for channel in range(3)
            ]
        ).T
        vertex_colors = np.repeat(pigment[:, None], len(cross), axis=1)
        if fancy and structured:
            if secondary_type == "H":
                flag = inverted[
                    np.clip(
                        (np.floor(x) if last and not first else np.ceil(x)).astype(int),
                        0,
                        len(ca) - 1,
                    )
                ]
                mask = np.where(
                    flag[:, None], (normals[:, 0] < -0.999)[None], backs[None]
                )
            else:
                mask = np.broadcast_to(sides, vertex_colors.shape[:2])
            original = vertex_colors[mask]
            value = original.max(axis=1, keepdims=True)
            saturation = (value - original.min(axis=1, keepdims=True)) / np.maximum(
                value, 1e-8
            )
            amount = np.maximum(saturation - 0.4, 0) / np.maximum(saturation, 1e-8)
            vertex_colors[mask] = value - (value - original) * amount
        owners = indices[np.clip(np.floor(x + 0.5).astype(int), 0, len(ca) - 1)]
        count = len(cross)
        a = (np.arange(len(x) - 1)[:, None] * count + np.arange(count)).ravel()
        b = (
            np.arange(len(x) - 1)[:, None] * count + (np.arange(count) + 1) % count
        ).ravel()
        faces = np.vstack((np.c_[a + count, a, b + count], np.c_[a, b, b + count]))
        if arrow:
            faces = faces[~np.any(faces < count, axis=1)]
            shoulders = np.vstack(
                (
                    np.c_[
                        np.arange(count),
                        np.arange(count) + count,
                        (np.arange(count) + 1) % count,
                    ],
                    np.c_[
                        (np.arange(count) + 1) % count,
                        np.arange(count) + count,
                        (np.arange(count) + 1) % count + count,
                    ],
                )
            )
            meshes.append(
                Mesh(
                    vertices[:2].reshape(-1, 3),
                    np.tile(-unit(tangent[0]), (2 * count, 1)),
                    np.tile(pigment[0], (2 * count, 1)),
                    shoulders,
                    np.repeat(owners[:2], count),
                )
            )
        meshes.append(
            Mesh(
                vertices.reshape(-1, 3),
                n.reshape(-1, 3),
                vertex_colors.reshape(-1, 3),
                faces,
                np.repeat(owners, count),
            )
        )
        for enabled, k, sign in (
            (j == 0 or (structured and first and kind != "ellipse"), 0, -1),
            (j == len(ca) - 1 or (structured and last and kind != "ellipse"), -1, 1),
        ):
            if not enabled:
                continue
            rim = vertices[k]
            cap_normal = sign * unit(tangent[k])
            if j == 0 and k == 0 and structured:
                rim = path[k] + e1[k] * cross[:, :1] + e2[k] * cross[:, 1:2]
            if not structured:
                t = np.linspace(0, 1, 6)
                radius = np.sqrt(1 - t * t)
                cap_vertices = (
                    path[k]
                    + radius[:, None, None] * (rim - path[k])
                    + coil * t[:, None, None] * cap_normal
                )
                cap_normals = unit(
                    radius[:, None, None] * n[k] + t[:, None, None] * cap_normal
                )
                a = (np.arange(5)[:, None] * count + np.arange(count)).ravel()
                b = (
                    np.arange(5)[:, None] * count + (np.arange(count) + 1) % count
                ).ravel()
                cap_faces = np.vstack(
                    (np.c_[a, b, a + count], np.c_[b, b + count, a + count])
                )
                meshes.append(
                    Mesh(
                        cap_vertices.reshape(-1, 3),
                        cap_normals.reshape(-1, 3),
                        np.tile(pigment[k], (6 * count, 1)),
                        cap_faces,
                        np.full(6 * count, owners[k]),
                    )
                )
                continue
            cap_faces = np.c_[
                np.full(count, count), np.arange(count), (np.arange(count) + 1) % count
            ]
            meshes.append(
                Mesh(
                    np.vstack((rim, path[k])),
                    np.tile(cap_normal, (count + 1, 1)),
                    np.tile(pigment[k], (count + 1, 1)),
                    cap_faces,
                    np.full(count + 1, owners[k]),
                )
            )
    for mesh in meshes:
        triangle = mesh.vertices[mesh.faces]
        face_normal = np.cross(
            triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0]
        )
        reverse = (
            np.sum(face_normal * mesh.normals[mesh.faces].mean(axis=1), axis=1) < 0
        )
        mesh.faces[reverse] = mesh.faces[reverse, ::-1]
    return meshes
