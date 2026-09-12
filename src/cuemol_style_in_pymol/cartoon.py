"""CueMol Ribbon2Renderer cylinder, sheet, and weighted coil splines."""

import itertools
from functools import lru_cache

import numpy as np
from scipy.interpolate import CubicSpline

from .mesh import unit


@lru_cache(maxsize=128)
def fit_operator(count, rho, weights):
    nodes = np.linspace(-1, 1, max(4, count))
    basis = CubicSpline(nodes, np.eye(len(nodes)), bc_type="natural")
    design = basis(np.linspace(-1, 1, count))
    weights = np.asarray(weights)
    weighted = design * weights[:, None]
    mid, half = (nodes[1:] + nodes[:-1]) / 2, np.diff(nodes) / 2
    sample = (mid[:, None] + half[:, None] * [-1, 1] / np.sqrt(3)).ravel()
    curvature = basis(sample, 2)
    penalty = curvature.T @ (np.repeat(half, 2)[:, None] * curvature)
    scale = np.sum(weighted**2, axis=0).max()
    strength = 10.0**rho * scale / np.diag(penalty).max()
    normal = weighted.T @ weighted + strength * penalty
    normal += np.eye(len(nodes)) * scale * (1 + strength) * 10 * np.finfo(float).eps
    return (nodes + 1) * (count - 1) / 2, np.linalg.solve(normal, weighted.T * weights)


def fit(points, rho, weights=None):
    weights = (1.0,) * len(points) if weights is None else tuple(weights)
    nodes, operator = fit_operator(len(points), rho, weights)
    return CubicSpline(nodes, operator @ points, bc_type="natural")


def build(segment, atoms, coords, colors, profile, axial, detail):
    from .geometry import frames, sweep

    indices = np.array([p for p, _ in segment])
    ca = coords[indices]
    secondary = np.array([atoms[i].ss for i in indices])
    boundaries = np.r_[0, np.flatnonzero(secondary[1:] != secondary[:-1]) + 1, len(ca)]
    runs = list(itertools.pairwise(boundaries))
    splines, meshes = {}, []
    rounded = profile.section == "ellipse"

    def samples(start, end):
        return np.linspace(start, end, max(1, int(np.floor((end - start) * axial))) + 1)

    def emit(path, tangent, hints, width, thickness, owners, kind):
        side, up = frames(path, hints, tangent)
        meshes.append(
            sweep(
                path,
                width,
                thickness,
                hints,
                colors[owners],
                owners,
                kind,
                detail,
                frame=(side, up),
            )
        )

    for start, end in runs:
        kind = secondary[start]
        if kind not in ("H", "S") or end - start < 2:
            continue
        first, last = max(0, start - 1), min(len(ca), end + 1)
        points = ca[first:last]
        spline = fit(points, 1.0 if kind == "S" and rounded else 3.0)
        splines[start] = (first, spline)
        begin, finish = (
            (start - 0.5 if first < start else start),
            (end - 0.5 if last > end else end - 1),
        )
        x = samples(begin, finish)
        path, tangent = spline(x - first), spline(x - first, 1)
        owners = indices[np.clip(np.floor(x + 0.5).astype(int), start, end - 1)]
        if kind == "H":
            radius = (
                np.linalg.norm(points - spline(np.arange(len(points))), axis=1).mean()
                + 0.2
            )
            hints = np.empty_like(path)
            previous = path[0] - ca[start]
            for j, velocity in enumerate(tangent):
                thin = unit(np.cross(previous, velocity))
                previous = unit(np.cross(velocity, thin))
                hints[j] = previous
            emit(path, tangent, hints, radius, radius, owners, "ellipse")
        else:
            normals = []
            normal_first = 2 if first < start else 1
            normal_last = len(points) - 2 if last > end else len(points) - 1
            for j in range(normal_first, normal_last + 1):
                residue = segment[first + j][1]
                if secondary[first + j] == "S" and "O" in residue and "C" in residue:
                    normal = unit(coords[residue["O"]] - coords[residue["C"]])
                elif 0 < j < len(points) - 1:
                    normal = unit(
                        np.cross(points[j] - points[j - 1], points[j + 1] - points[j])
                    )
                else:
                    normal = np.array([1.0, 0.0, 0.0])
                if normals and np.dot(normal, normals[-1]) < 0:
                    normal = -normal
                normals.append(normal)
            if len(normals) > 1 and normal_last < len(points) - 1:
                hints = unit(fit(np.asarray(normals), 5.0)(x - first - normal_first))
            else:
                hints = np.tile(
                    unit(np.sum(normals, axis=0)) if normals else [1, 0, 0], (len(x), 1)
                )
            width = np.full(len(x), 1.4)
            arrow = x >= finish - 1
            width[arrow] = 0.2 + (1.4 * 1.8 - 0.2) * (finish - x[arrow])
            emit(path, tangent, hints, width, 0.2, owners, profile.section)

    for number, (start, end) in enumerate(runs):
        if secondary[start] in ("H", "S") and end - start >= 2:
            continue
        points, weights = list(ca[start:end]), [1.0] * (end - start)
        render_start, render_end = 0, len(points) - 1
        if start > 0:
            previous = runs[number - 1][0]
            if secondary[start - 1] == "S" and previous in splines:
                first, spline = splines[previous]
                points[:0] = [spline(start - 1.5 - first), spline(start - 0.5 - first)]
                weights[:0] = [10.0, 10.0]
                render_start, render_end = 1, render_end + 2
            else:
                points.insert(0, ca[start - 1])
                weights.insert(0, 10.0)
                weights[1] = 10.0
                render_end += 1
        if end < len(ca):
            if secondary[end] == "S" and end in splines:
                first, spline = splines[end]
                points.extend([spline(end - 0.5 - first), spline(end + 0.5 - first)])
                weights.extend([10.0, 10.0])
                render_end += 1
            else:
                weights[-1] = 10.0
                points.append(ca[end])
                weights.append(10.0)
                render_end += 1
        else:
            weights[-1] = 10.0
        if len(points) < 2:
            continue
        spline = fit(np.asarray(points), -2.0 if rounded else -1.0, weights)
        t = samples(render_start, render_end)
        path, tangent = spline(t), spline(t, 1)
        owners = indices[
            np.clip(
                np.floor(t - render_start + start + 0.5).astype(int), start, end - 1
            )
        ]
        emit(
            path,
            tangent,
            np.tile([1.0, 0.0, 0.0], (len(t), 1)),
            0.2,
            0.2,
            owners,
            "ellipse",
        )
    return meshes
