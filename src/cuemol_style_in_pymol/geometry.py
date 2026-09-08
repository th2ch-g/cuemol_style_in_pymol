"""Molecular mesh construction inspired by CueMol's documented sections.

The implementation uses Catmull-Rom interpolation and transported peptide
frames, not CueMol's smoothing-spline implementation. Widths follow the
public DefaultRibbon/Fancy1Ribbon style definitions (angstrom units).
"""

import colorsys
from functools import lru_cache
import xml.etree.ElementTree as ET

import numpy as np

from .mesh import Mesh, merge, unit
from .presets import (
    COIL_COLOR,
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
    t = np.linspace(0, 2 * np.pi, detail, endpoint=False)
    points = np.c_[np.cos(t), np.sin(t)]
    if kind == "fancy":
        # A thin middle with rounded rails at both ribbon margins.
        points[:, 1] *= 0.3 + 0.7 * np.abs(points[:, 0]) ** 6
    tangents = np.roll(points, -1, axis=0) - np.roll(points, 1, axis=0)
    return points, unit(np.c_[tangents[:, 1], -tangents[:, 0]])


def interpolate(points, samples):
    p = np.asarray(points, float)
    if len(p) == 1:
        return p.copy(), np.array([0.0])
    x = np.arange((len(p) - 1) * samples + 1) / samples
    i = np.minimum(x.astype(int), len(p) - 2)
    t = (x - i)[:, None]
    a, b, c, d = [p[np.clip(i + j, 0, len(p) - 1)] for j in (-1, 0, 1, 2)]
    return 0.5 * (
        (2 * b)
        + (-a + c) * t
        + (2 * a - 5 * b + 4 * c - d) * t * t
        + (-a + 3 * b - 3 * c + d) * t * t * t
    ), x


def frames(path, hints):
    tangent = unit(np.gradient(path, axis=0))
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
    path, widths, thickness, hints, colors, owners, kind, detail, back=False, front=None
):
    if len(path) < 2:
        return Mesh([], [], [], [], [])
    side, up = frames(path, hints)
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
    if back:
        # Frame signs can depend on preceding strands and loops. Color helix
        # undersides by the physical inward direction without twisting the mesh.
        polarity = (
            np.where(np.sum(up * front, axis=1) < 0, -1, 1)
            if front is not None
            else np.ones(len(path))
        )
        mask = polarity[:, None] * shape[None, :, 1] < -0.1
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
        normal = unit(path[1] - path[0]) if idx == 0 else unit(path[-1] - path[-2])
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
    return sweep(
        np.array([a, (a + b) / 2, b]),
        radius,
        radius,
        np.tile(h, (3, 1)),
        [colors[0], colors[0], colors[1]],
        [owners[0], owners[0], owners[1]],
        "ellipse",
        detail,
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
        pivot = r.get("CA", r.get("P", r.get("C4'", r.get("C4*"))))
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
        hints = []
        for j, (p, r) in enumerate(segment):
            if "C" in r and "O" in r:
                hints.append(coords[r["O"]] - coords[r["C"]])
            else:
                d = ca[min(j + 1, len(ca) - 1)] - ca[max(j - 1, 0)]
                hints.append(np.eye(3)[np.argmin(np.abs(d))])
        path, x = interpolate(ca, axial)
        # Align peptide frames before interpolation: beta-strand carbonyls
        # alternate sides and otherwise make the interpolated frame collapse.
        aligned, _ = frames(ca, np.asarray(hints))
        hi, _ = interpolate(aligned, axial)
        pick = np.clip(np.floor(x + 0.5).astype(int), 0, len(indices) - 1)
        owner = indices[pick]
        ss = np.array([atoms[i].ss for i in owner])
        seq = [atoms[i].ss for i in indices]
        front = None
        if profile.back and representation == "ribbon":
            front, _ = interpolate(helix_outward(ca, np.asarray(hints), seq), axial)
            front[ss != "H"] = 0
        width = np.full(len(path), 0.25)
        thickness = np.full(len(path), 0.25)
        if representation not in ("tube", "nucleic"):
            width[ss == "H"] = 1.3 if profile.section == "fancy" else 1.2
            width[ss == "S"] = 1.2
            thickness[(ss == "H") | (ss == "S")] = 0.2
            # Taper the final residue of each sheet toward the C terminus.
            for j, s in enumerate(seq):
                if s == "S" and (j == len(seq) - 1 or seq[j + 1] != "S"):
                    tip = j if j == len(seq) - 1 else j + 0.5
                    region = (x >= max(0, tip - 1.2)) & (x <= tip)
                    width[region] = np.maximum(0.025, 1.92 * (tip - x[region]) / 1.2)
        if representation == "cartoon":
            start = 0
            while start < len(seq):
                end = start + 1
                while end < len(seq) and seq[end] == seq[start]:
                    end += 1
                if seq[start] == "H" and end - start >= 3:
                    part = ca[start:end]
                    center = part.mean(axis=0)
                    axis = np.linalg.svd(part - center, full_matrices=False)[2][0]
                    if np.dot(axis, part[-1] - part[0]) < 0:
                        axis = -axis
                    projected = center + np.outer((part - center) @ axis, axis)
                    mask = (x >= start) & (x <= end - 1)
                    for k in range(3):
                        path[mask, k] = np.interp(
                            x[mask], np.arange(start, end), projected[:, k]
                        )
                    width[mask] = thickness[mask] = 1.6
                start = end
        # Sections change only at secondary-structure boundaries.
        kinds = np.where(np.isin(ss, ["H", "S"]), profile.section, "ellipse")
        if profile.section == "fancy":
            kinds[ss == "S"] = "rectangle"
        if representation in ("tube", "nucleic"):
            kinds[:] = "ellipse"
        if representation == "cartoon":
            kinds[ss == "H"] = "ellipse"
        begin = 0
        while begin < len(path) - 1:
            end = begin + 1
            while end < len(path) - 1 and kinds[end] == kinds[begin]:
                end += 1
            sl = slice(begin, end + 1)
            meshes.append(
                sweep(
                    path[sl],
                    width[sl],
                    thickness[sl],
                    hi[sl],
                    colors[owner[sl]],
                    owner[sl],
                    str(kinds[begin]),
                    detail,
                    profile.back and kinds[begin] != "ellipse",
                    front[sl] if front is not None else None,
                )
            )
            begin = end
    return merge(meshes)


def nucleic_bases(atoms, coords, colors, detail):
    residues = {}
    for i, a in enumerate(atoms):
        if a.kind == "nucleic":
            residues.setdefault((a.segi, a.chain, a.resi), []).append(i)
    meshes = []
    from scipy.spatial import ConvexHull

    for indices in residues.values():
        base = [
            i
            for i in indices
            if "'" not in atoms[i].name
            and "*" not in atoms[i].name
            and atoms[i].element in ("C", "N")
            and atoms[i].name not in ("C5M",)
        ]
        if len(base) < 3:
            continue
        p = coords[base]
        center = p.mean(axis=0)
        axes = np.linalg.svd(p - center, full_matrices=False)[2]
        flat = (p - center) @ axes[:2].T
        try:
            ring = ConvexHull(flat).vertices
        except Exception:
            continue
        ringp = center + flat[ring] @ axes[:2]
        owner = base[0]
        n = axes[2]
        count = len(ring)
        v = np.vstack((ringp + 0.15 * n, ringp - 0.15 * n))
        f = []
        for j in range(1, count - 1):
            f.extend(((0, j, j + 1), (count, count + j + 1, count + j)))
        for j in range(count):
            k = (j + 1) % count
            f.extend(((j, k, count + j), (k, count + k, count + j)))
        faces = np.asarray(f)
        normals = np.cross(
            v[faces[:, 1]] - v[faces[:, 0]], v[faces[:, 2]] - v[faces[:, 0]]
        )
        reverse = np.sum(normals * (v[faces].mean(axis=1) - center), axis=1) < 0
        faces[reverse] = faces[reverse, ::-1]
        normals[reverse] *= -1
        # Duplicate triangle vertices for flat base faces and sharp slab rims.
        v = v[faces].reshape(-1, 3)
        normals = np.repeat(unit(normals), 3, axis=0)
        faces = np.arange(len(v)).reshape(-1, 3)
        meshes.append(
            Mesh(
                v,
                normals,
                np.tile(colors[owner], (len(v), 1)),
                faces,
                np.full(len(v), owner),
            )
        )
        anchor = next(
            (i for i in indices if atoms[i].name in ("C1'", "C1*", "P")), None
        )
        if anchor is not None:
            meshes.append(
                bond(
                    coords[anchor],
                    center,
                    0.15,
                    [colors[owner]] * 2,
                    [owner] * 2,
                    detail,
                )
            )
    return merge(meshes)


def surface_mesh(model, colors, coords, quality):
    """Use PyMOL's SES algorithm in an isolated, headless instance."""
    import pymol2
    from scipy.spatial import cKDTree

    with pymol2.PyMOL() as p:
        c = p.cmd
        c.load_model(model, "surface_source")
        c.hide("everything")
        c.show("surface")
        c.set("surface_quality", {"low": 0, "medium": 1, "high": 2}[quality])
        c.set("surface_solvent", 0)
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
            radius = a.vdw
        elif representation == "sticks":
            radius = 0.18
        else:
            radius = 0.12 if a.element == "H" else 0.3
        parts.append(sphere(coords[i], radius, colors[i], i, detail))
    if representation != "cpk":
        for i, j in bonds:
            if i in chosen and j in chosen:
                parts.append(
                    bond(
                        coords[i],
                        coords[j],
                        0.14,
                        [colors[i], colors[j]],
                        [i, j],
                        detail,
                    )
                )
    return merge(parts)
