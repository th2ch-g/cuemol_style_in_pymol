"""Image-space contours from visible surface continuity, independent of mesh edges."""

from collections import Counter, defaultdict
from itertools import combinations, pairwise

import numpy as np

from ._edtsurf import sample_image, screen_cracks, trace_paths
from .mesh import merge, unit


def surface_image(mesh, matrices, budget, samples=3):
    from .sampling import clip_mesh, project

    w, h = (matrices[2][2:] * samples).astype(int)
    if w * h * 40 > budget:
        raise ValueError("Contour surface samples exceed cache_mb")
    clipped = clip_mesh(mesh, matrices)
    xy, ndc, eye = project(clipped.vertices, matrices)
    inverse_w = 1 / np.maximum(matrices[1][3, 2] * eye[:, 2] + matrices[1][3, 3], 1e-8)
    attributes = np.c_[
        -eye[:, 2], clipped.normals @ matrices[0][:3, :3].T, clipped.colors
    ]
    return sample_image(
        np.c_[xy * samples, ndc[:, 2], inverse_w],
        clipped.faces,
        attributes,
        int(w),
        int(h),
    )


def prune(begin, end, flags):
    active = np.ones(len(flags), bool)
    for _ in range(8):
        chains = trace_paths(begin, end, active)
        kept, interior = [], []
        for points, ids, d0, d1 in chains:
            values = flags[ids]
            count = np.count_nonzero(values & 32)
            kept.append(
                bool(
                    np.any((values & 7) == 1)
                    or count >= 6
                    or (d0 >= 3 and d1 >= 3 and count)
                )
            )
            interior.append(bool(count))
        while True:
            support, inside = set(), set()
            for (points, _, _, _), keep, inner in zip(chains, kept, interior):
                if keep:
                    support.update((tuple(points[0]), tuple(points[-1])))
                    if inner:
                        inside.update((tuple(points[0]), tuple(points[-1])))
            changed = False
            for i, (points, _, _, _) in enumerate(chains):
                ends = {tuple(points[0]), tuple(points[-1])}
                if not kept[i] and len(ends) == 2 and ends <= support and ends & inside:
                    kept[i] = interior[i] = True
                    changed = True
            if not changed:
                break
        previous = active.copy()
        endpoints = Counter(
            tuple(p)
            for (points, _, _, _), keep in zip(chains, kept)
            if keep and not np.array_equal(points[0], points[-1])
            for p in (points[0], points[-1])
        )
        for (points, ids, _, _), keep in zip(chains, kept):
            if not keep:
                active[ids] = False
                continue
            if np.array_equal(points[0], points[-1]):
                continue
            values = flags[ids]
            strong = np.count_nonzero(values & 32) >= 6
            for reverse in (False, True):
                if endpoints[tuple(points[-1 if reverse else 0])] > 1:
                    continue
                order = ids[::-1] if reverse else ids
                remove = []
                for i in order:
                    weak = (flags[i] & 7) == 3 and not flags[i] & 32
                    if not weak or (strong and not flags[i] & 64):
                        break
                    remove.append(i)
                if len(remove) != len(ids):
                    active[remove] = False
        if np.array_equal(previous, active):
            return chains
    return trace_paths(begin, end, active)


def simplify(points, epsilon=1.2):
    closed = np.array_equal(points[0, :2], points[-1, :2])
    keep = np.zeros(len(points), bool)
    keep[[0, -1]] = True
    stack = [(0, len(points) - 1)]
    if closed and len(points) >= 4:
        unique = points[:-1]
        a = np.argmax(np.sum((unique[:, :2] - unique[0, :2]) ** 2, axis=1))
        b = np.argmax(np.sum((unique[:, :2] - unique[a, :2]) ** 2, axis=1))
        lo, hi = sorted((a, b))
        unique = np.roll(unique, -lo, axis=0)
        points = np.vstack((unique, unique[0]))
        keep[hi - lo] = True
        stack = [(0, hi - lo), (hi - lo, len(points) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        delta = points[b, :2] - points[a, :2]
        length = np.linalg.norm(delta)
        values = points[a + 1 : b, :2] - points[a, :2]
        distances = (
            np.abs(values[:, 0] * delta[1] - values[:, 1] * delta[0]) / length
            if length > 0
            else np.linalg.norm(values, axis=1)
        )
        index = int(np.argmax(distances))
        if distances[index] > epsilon:
            pivot = a + index + 1
            keep[pivot] = True
            stack.extend(((a, pivot), (pivot, b)))
    return points[keep]


def rewire(chains, width):
    """Reconnect a contour when pixel continuity leaves a short junction gap."""
    from scipy.spatial import cKDTree

    chains = list(chains)
    radius = min(8, max(2, int(width / 2 + 0.5) + 2))
    for _ in range(4):
        mids = (
            np.concatenate([(p[:-1] + p[1:]) / 2 for p, _, _, _ in chains])
            if chains
            else np.empty((0, 2))
        )
        source = (
            np.concatenate(
                [np.full(len(ids), i) for i, (_, ids, _, _) in enumerate(chains)]
            )
            if chains
            else np.empty(0, int)
        )
        edge = (
            np.concatenate([np.arange(len(ids)) for _, ids, _, _ in chains])
            if chains
            else np.empty(0, int)
        )
        tree = cKDTree(mids)
        candidates = []
        for i, (points, ids, d0, d1) in enumerate(chains):
            if len(points) < 6 or np.array_equal(points[0], points[-1]):
                continue
            for side, degree in enumerate((d0, d1)):
                if degree > 1:
                    continue
                ordered = points[::-1] if side else points
                origin = ordered[0]
                outward = unit(origin - ordered[min(8, len(ids))])
                nearby = np.asarray(tree.query_ball_point(origin, radius + 0.5), int)
                delta = mids[nearby] - origin
                distance = np.linalg.norm(delta, axis=1)
                valid = (
                    (source[nearby] != i)
                    & (distance > 0.5)
                    & (delta @ outward >= 0.7071 * distance)
                )
                if not valid.any():
                    continue
                winner = nearby[np.argmin(np.where(valid, distance, np.inf))]
                target, hit = source[winner], edge[winner]
                path = chains[target][0]
                if np.array_equal(path[0], path[-1]):
                    continue
                indices = np.arange(
                    max(8, hit - 2 * radius), min(len(path) - 8, hit + 2 * radius + 1)
                )
                if not len(indices):
                    continue
                before = unit(path[indices - 8] - path[indices])
                after = unit(path[indices + 8] - path[indices])
                through = -np.sum(before * after, axis=1)
                chosen = np.argmin(through)
                high = after[chosen] @ outward >= before[chosen] @ outward
                continuation = (after if high else before)[chosen] @ outward
                if (
                    continuation >= 0.82
                    and through[chosen] < 0.82
                    and continuation >= through[chosen] + 0.1
                ):
                    candidates.append((i, side, target, indices[chosen], high))
        if not candidates:
            break
        touched = set()
        for i, side, target, j, high in candidates:
            if i in touched or target in touched:
                continue
            touched.update((i, target))
            p, ids, d0, d1 = chains[target]
            low = (p[: j + 1], ids[:j], d0, 3)
            upper = (p[j:], ids[j:], 3, d1)
            q, qids, qd0, qd1 = chains[i]
            if not side:
                q, qids, qd0 = q[::-1], qids[::-1], qd1
            r, rids, rd0, rd1 = upper if high else low
            if not high:
                r, rids, rd1 = r[::-1], rids[::-1], rd0
            chains[i] = (np.vstack((q, r)), np.r_[qids, qids[-1], rids], qd0, rd1)
            chains[target] = low if high else upper
    return chains


def weave(chains, flags, depths, pixel):
    """Pair junction ends once, then assemble each continuous contour once."""
    from scipy.spatial import cKDTree

    ends, positions = [], []
    for i, (points, ids, _, _) in enumerate(chains):
        if np.array_equal(points[0], points[-1]):
            continue
        for side in (0, 1):
            ordered = points[::-1] if side else points
            ends.append(
                (
                    2 * i + side,
                    unit(ordered[min(8, len(ids))] - ordered[0]),
                    ids[-1 if side else 0],
                )
            )
            positions.append(ordered[0])
    if not ends:
        return chains
    parent = list(range(len(ends)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in sorted(cKDTree(positions).query_pairs(2, p=np.inf)):
        parent[root(b)] = root(a)
    clusters = defaultdict(list)
    for i, end in enumerate(ends):
        clusters[root(i)].append(end)
    partner = {}
    for members in clusters.values():
        candidates = []
        for (a, da, ia), (b, db, ib) in combinations(members, 2):
            cosine = da @ db
            if cosine <= -0.82 and abs(depths[ia] - depths[ib]) <= 300 * pixel:
                candidates.append(((flags[ia] & 7) != (flags[ib] & 7), cosine, a, b))
        for _, _, a, b in sorted(candidates):
            if a not in partner and b not in partner:
                partner[a], partner[b] = b, a
    used, result = set(), []
    starts = [e for e, _, _ in ends if e not in partner] + [
        2 * i for i in range(len(chains))
    ]
    for start in starts:
        current, side = divmod(start, 2)
        if current in used:
            continue
        parts, edge_parts = [], []
        degree_start = chains[current][2 + side]
        closed = False
        while current not in used:
            used.add(current)
            points, ids, d0, d1 = chains[current]
            if side:
                points, ids, d0, d1 = points[::-1], ids[::-1], d1, d0
            if parts:
                gap = not np.array_equal(parts[-1][-1], points[0])
                if gap:
                    edge_parts.append(edge_parts[-1][-1:])
                parts.append(points if gap else points[1:])
            else:
                parts.append(points)
            edge_parts.append(ids)
            degree_end = d1
            next_end = partner.get(2 * current + 1 - side)
            if next_end is None:
                break
            current, side = divmod(next_end, 2)
            if current in used:
                closed = True
        points, ids = np.vstack(parts), np.concatenate(edge_parts)
        if closed:
            if not np.array_equal(points[0], points[-1]):
                points = np.vstack((points, points[0]))
                ids = np.r_[ids, ids[-1]]
            degree_start = degree_end = 2
        result.append((points, ids, degree_start, degree_end))
    return result


def smooth(points, width):
    closed = np.array_equal(points[0, :2], points[-1, :2])
    delta = np.diff(points[:, :2], axis=0)
    cross = delta[:-1, 0] * delta[1:, 1] - delta[:-1, 1] * delta[1:, 0]
    keep = np.r_[
        True, (cross != 0) | (np.sum(delta[:-1] * delta[1:], axis=1) <= 0), True
    ]
    points = points[keep]
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(points[:, :2], axis=0), axis=1))]
    pins = np.zeros(len(points), bool)
    a = np.maximum(0, np.searchsorted(arc, arc[1:-1] - 4, side="right") - 1)
    b = np.minimum(len(points) - 1, np.searchsorted(arc, arc[1:-1] + 4))
    before = unit(points[1:-1, :2] - points[a, :2])
    after = unit(points[b, :2] - points[1:-1, :2])
    pins[1:-1] = np.sum(before * after, axis=1) < 0.5
    delta = np.diff(points, axis=0)
    counts = np.maximum(
        1, np.ceil(np.linalg.norm(delta[:, :2], axis=1) / max(4, width))
    ).astype(int)
    offsets = np.r_[0, np.cumsum(counts)]
    source = np.repeat(np.arange(len(delta)), counts)
    fraction = (np.arange(offsets[-1]) - np.repeat(offsets[:-1], counts)) / counts[
        source
    ]
    subdivided = points[source] + fraction[:, None] * delta[source]
    markers = np.zeros(offsets[-1] + 1, bool)
    markers[offsets] = pins
    points = np.vstack((subdivided, points[-1]))
    pins = markers

    def chaikin(part, loop):
        for _ in range(2):
            if len(part) < 3:
                break
            q = 0.75 * part[:-1] + 0.25 * part[1:]
            r = 0.25 * part[:-1] + 0.75 * part[1:]
            result = np.stack((q, r), axis=1).reshape(-1, part.shape[1])
            part = (
                np.vstack((result, result[0]))
                if loop
                else np.vstack((part[0], result, part[-1]))
            )
        return part

    if not pins.any():
        return simplify(chaikin(points, closed))
    if closed:
        pivot = np.flatnonzero(pins)[0]
        points = np.roll(points[:-1], -pivot, axis=0)
        points = np.vstack((points, points[0]))
        pins = np.roll(pins[:-1], -pivot)
        pins = np.r_[pins, pins[0]]
    cuts = np.unique(np.r_[0, np.flatnonzero(pins), len(points) - 1])
    parts = [chaikin(points[a : b + 1], False) for a, b in pairwise(cuts)]
    return simplify(np.vstack([parts[0]] + [part[1:] for part in parts[1:]]))


def contour_image(
    meshes,
    matrices,
    width,
    budget,
    outer_only=False,
    far_limit=np.inf,
    analytic=False,
    camera_depth=None,
):
    mesh = merge(meshes)
    surface, depth = surface_image(mesh, matrices, budget)
    begin, end, flags, z, owners = screen_cracks(
        surface,
        matrices[1],
        mesh.vertices @ matrices[0][:3, :3].T + matrices[0][:3, 3],
        mesh.faces,
        outer_only,
        far_limit,
        analytic,
    )
    chains = prune(begin, end, flags)
    pixel_width = max(3, width * matrices[1][1, 1] * surface.shape[0] / 2)
    if abs(matrices[1][3, 3]) < 0.5:
        center_depth = max(camera_depth or -matrices[0][2, 3], 1e-6)
        pixel_width = max(3, pixel_width / center_depth)
    pixel = 2 / (matrices[1][1, 1] * surface.shape[0])
    if abs(matrices[1][3, 3]) < 0.5:
        pixel *= center_depth
    chains = weave(rewire(chains, pixel_width), flags, z, pixel)
    polylines = []
    for points, ids, d0, d1 in chains:
        closed = np.array_equal(points[0], points[-1])
        if len(ids) < 12 and (closed or d0 < 3 or d1 < 3):
            continue
        direction = np.diff(points, axis=0)
        offset = owners[ids] - (points[:-1] + points[1:]) / 2
        outside = -np.sign(
            np.sum(
                np.sign(direction[:, 0] * offset[:, 1] - direction[:, 1] * offset[:, 0])
            )
        )
        vz = np.r_[z[ids[0]], (z[ids[:-1]] + z[ids[1:]]) / 2, z[ids[-1]]]
        if closed:
            vz[[0, -1]] = (z[ids[0]] + z[ids[-1]]) / 2
        points = smooth(np.c_[points, vz], pixel_width)
        polylines.append((points, outside, d0 >= 3, d1 >= 3))
    return surface, depth, polylines, pixel_width


def stroke_image(polylines, width, shape):
    """Rasterize round, joined bands with a half-pixel pad on the near side."""
    vertices, triangles = [], []

    def polygon(points):
        first = len(vertices)
        vertices.extend(points)
        triangles.extend(
            (first, first + i, first + i + 1) for i in range(1, len(points) - 1)
        )

    for points, outside, taper_start, taper_end in sorted(
        polylines, key=lambda p: -p[0][:, 2].min()
    ):
        if len(points) < 2:
            continue
        closed = np.array_equal(points[0, :2], points[-1, :2])
        lengths = np.linalg.norm(np.diff(points[:, :2], axis=0), axis=1)
        arc = np.r_[0, np.cumsum(lengths)]
        # Preserve the authored corners while filling long interpolation spans.
        positions = np.unique(
            np.concatenate(
                [arc, *[a + np.arange(6, b - a, 6) for a, b in pairwise(arc)]]
            )
        )
        points = np.c_[[np.interp(positions, arc, points[:, k]) for k in range(3)]].T
        direction = unit(np.diff(points[:, :2], axis=0))
        normals = np.c_[-direction[:, 1], direction[:, 0]]
        starts, ends = normals.copy(), normals.copy()
        for i in range(1, len(points) - 1):
            if direction[i - 1] @ direction[i] >= 0.9848:
                starts[i] = ends[i - 1] = unit(normals[i - 1] + normals[i])
        if closed and direction[-1] @ direction[0] >= 0.9848:
            starts[0] = ends[-1] = unit(normals[-1] + normals[0])
        pad, outer = min(width / 2, 1.5), max(width / 2, width - 1.5)
        left = np.full(
            len(points), outer if outside > 0 else pad if outside < 0 else width / 2
        )
        right = np.full(
            len(points), pad if outside > 0 else outer if outside < 0 else width / 2
        )
        distance = np.full(len(points), width)
        if taper_start:
            distance = np.minimum(distance, positions)
        if taper_end:
            distance = np.minimum(distance, positions[-1] - positions)
        k = np.clip(distance / width, 0, 1)
        k = k * k * (3 - 2 * k)
        left = width / 2 + (left - width / 2) * k
        right = width - left
        for i, normal in enumerate(normals):
            a, b = points[i], points[i + 1]
            polygon(
                np.array(
                    [
                        [*(a[:2] + starts[i] * left[i]), a[2]],
                        [*(a[:2] - starts[i] * right[i]), a[2]],
                        [*(b[:2] - ends[i] * right[i + 1]), b[2]],
                        [*(b[:2] + ends[i] * left[i + 1]), b[2]],
                    ]
                )
            )
        for i in range(len(points)):
            if closed and i == len(points) - 1:
                continue
            endpoint = i in (0, len(points) - 1)
            if endpoint and not closed:
                if (i == 0 and taper_start) or (i and taper_end):
                    continue
                normal = normals[0 if i == 0 else -1]
                center = points[i, :2] + normal * (left[i] - right[i]) / 2
                start = np.arctan2(normal[1], normal[0])
                span = np.pi if i == 0 else -np.pi
                radius = width / 2
            else:
                previous = normals[i - 1]
                following = normals[i]
                if direction[i - 1] @ direction[i] >= 0.9848:
                    continue
                cross = previous[0] * following[1] - previous[1] * following[0]
                if abs(cross) < 1e-8:
                    continue
                sign = -1 if cross > 0 else 1
                start = np.arctan2(previous[1] * sign, previous[0] * sign)
                finish = np.arctan2(following[1] * sign, following[0] * sign)
                span = (finish - start + np.pi) % (2 * np.pi) - np.pi
                center = points[i, :2]
                radius = right[i] if sign < 0 else left[i]
            count = max(2, int(np.ceil(abs(span) * radius)))
            theta = start + np.linspace(0, span, count + 1)
            xy = center + radius * np.c_[np.cos(theta), np.sin(theta)]
            polygon(np.c_[np.vstack((center, xy)), np.full(count + 2, points[i, 2])])
    if not vertices:
        return np.zeros(shape, np.float32)
    vertices = np.asarray(vertices)
    # Later bands win at crossings, matching far-to-near stroke painting.
    priority = -np.arange(len(vertices), dtype=float) / max(len(vertices), 1)
    image, _ = sample_image(
        np.c_[vertices[:, :2], priority, np.ones(len(vertices))],
        np.asarray(triangles),
        vertices[:, 2:3],
        shape[1],
        shape[0],
    )
    return image[:, :, 0]
