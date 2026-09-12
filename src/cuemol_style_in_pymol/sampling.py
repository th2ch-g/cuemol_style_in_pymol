"""Camera-dependent native material samples, prepared outside GL callbacks."""

import numpy as np

from ._edtsurf import pencil, visibility, visible_samples
from .materials import bake, drawing_tone, native_blend
from .mesh import Mesh, unit


def clipping(view):
    """Match PyMOL's safe near plane and minimum one-angstrom slab."""
    near, far = float(view[15]), float(view[16])
    if far - near < 1:
        near, far = (near + far - 1) / 2, (near + far + 1) / 2
    if near < 1:
        near, far = 1, max(far, 2)
    return near, far


def camera(cmd, width=0, height=0, *, ray=False):
    from .export import view_matrix

    view = np.asarray(cmd.get_view())
    vw, vh = cmd.get_viewport()
    if width and not height:
        height = round(width * vh / max(vw, 1))
    elif height and not width:
        width = round(height * vw / max(vh, 1))
    width, height = int(width or vw), int(height or vh)
    if width <= 0 or height <= 0:
        raise ValueError("A positive viewport is required for material sampling")
    tangent = np.tan(np.deg2rad(cmd.get_setting_float("field_of_view")) / 2)
    near, far = clipping(view)
    projection = np.zeros((4, 4))
    orthoscopic = cmd.get_setting_int("orthoscopic")
    if ray and cmd.get_setting_int("ray_orthoscopic") >= 0:
        orthoscopic = cmd.get_setting_int("ray_orthoscopic")
    if orthoscopic:
        scale = 1 / max(-view[11] * tangent, 1e-8)
        projection[0, 0], projection[1, 1] = scale * height / width, scale
        projection[2, 2], projection[2, 3] = (
            -2 / (far - near),
            -(far + near) / (far - near),
        )
        projection[3, 3] = 1
    else:
        # PyMOL 3.1 passes 2*tan(fov/2), rather than fov, to glm::perspective.
        if not ray:
            tangent = np.tan(tangent)
        projection[0, 0], projection[1, 1] = height / width / tangent, 1 / tangent
        projection[2, 2], projection[2, 3] = (
            -(far + near) / (far - near),
            -2 * far * near / (far - near),
        )
        projection[3, 2] = -1
    return view_matrix(cmd), projection, np.array([0, 0, width, height])


def project(vertices, matrices):
    modelview, projection, viewport = matrices
    eye = vertices @ modelview[:3, :3].T + modelview[:3, 3]
    clip = np.c_[eye, np.ones(len(eye))] @ projection.T
    divisor = np.maximum(clip[:, 3:4], 1e-8)
    ndc = clip[:, :3] / divisor
    return (
        np.c_[(ndc[:, 0] + 1) * viewport[2] / 2, (1 - ndc[:, 1]) * viewport[3] / 2],
        ndc,
        eye,
    )


def cancel_native_fog(colors, mesh, matrices, background, fog_range):
    if fog_range is None:
        return colors
    eye = project(mesh.vertices, matrices)[2]
    factor = np.clip(
        (fog_range[1] + eye[:, 2]) / max(fog_range[1] - fog_range[0], 1e-8), 1e-6, 1
    )
    return np.clip(background + (colors - background) / factor[:, None], 0, 1)


def clip_mesh(mesh, matrices):
    transform = matrices[1] @ matrices[0]
    vertices, normals, colors, owners = (
        mesh.vertices,
        mesh.normals,
        mesh.colors,
        mesh.owners,
    )
    faces = mesh.faces
    for axis in range(3):
        for sign in (-1, 1):
            homogeneous = np.c_[vertices, np.ones(len(vertices))] @ transform.T
            distance = homogeneous[:, 3] + sign * homogeneous[:, axis]
            inside = distance[faces] >= 0
            crossing = inside.any(axis=1) & ~inside.all(axis=1)
            additions = []
            triangles = []
            for face in faces[crossing]:
                polygon = []
                for a, b in zip(face, np.roll(face, -1)):
                    if distance[a] >= 0:
                        polygon.append(int(a))
                    if (distance[a] >= 0) != (distance[b] >= 0):
                        weight = distance[a] / (distance[a] - distance[b])
                        additions.append(
                            (
                                vertices[a] + weight * (vertices[b] - vertices[a]),
                                normals[a] + weight * (normals[b] - normals[a]),
                                colors[a] + weight * (colors[b] - colors[a]),
                                owners[a] if weight < 0.5 else owners[b],
                            )
                        )
                        polygon.append(len(vertices) + len(additions) - 1)
                triangles.extend(
                    (polygon[0], polygon[i], polygon[i + 1])
                    for i in range(1, len(polygon) - 1)
                )
            faces = faces[inside.all(axis=1)]
            if additions:
                positions, ns, cs, ids = zip(*additions)
                vertices = np.concatenate((vertices, positions))
                normals = np.concatenate((normals, ns))
                colors = np.concatenate((colors, cs))
                owners = np.concatenate((owners, ids))
                faces = np.vstack((faces, triangles)).astype(np.uint32)
    return Mesh(vertices, normals, colors, faces, owners, mesh.opacity)


def subdivide(mesh, matrices, budget, max_pixels=1 / 3):
    """Bisect projected longest edges until every visible edge meets the bound."""
    mesh = clip_mesh(mesh, matrices)
    vertices, normals, colors, owners = (
        mesh.vertices,
        mesh.normals,
        mesh.colors,
        mesh.owners,
    )
    faces = mesh.faces.copy()
    finished = []
    count = 0
    for _ in range(40):
        pixels, ndc, _ = project(vertices, matrices)
        triangles = ndc[faces]
        visible = np.all(triangles.max(axis=1) >= -1, axis=1) & np.all(
            triangles.min(axis=1) <= 1, axis=1
        )
        faces = faces[visible]
        if not len(faces):
            break
        points = pixels[faces]
        lengths = np.sum((points - np.roll(points, -1, axis=1)) ** 2, axis=2)
        edge = lengths.argmax(axis=1)
        split = lengths.max(axis=1) > max_pixels**2
        finished.append(faces[~split])
        count += int((~split).sum())
        if not split.any():
            faces = faces[:0]
            break
        faces, edge = faces[split], edge[split]
        projected_bytes = (count + len(faces) * 2) * 144 + (
            len(vertices) + len(faces)
        ) * 40
        if projected_bytes > budget:
            raise ValueError(
                "Camera-dependent material samples exceed cache_mb; increase cache_mb or reduce the image size"
            )
        rows = np.arange(len(faces))
        a, b, c = (
            faces[rows, edge],
            faces[rows, (edge + 1) % 3],
            faces[rows, (edge + 2) % 3],
        )
        low, high = (
            np.minimum(a, b).astype(np.uint64),
            np.maximum(a, b).astype(np.uint64),
        )
        codes, inverse = np.unique((low << np.uint64(32)) | high, return_inverse=True)
        endpoints = np.c_[codes >> np.uint64(32), codes & np.uint64(0xFFFFFFFF)].astype(
            int
        )
        mid = (len(vertices) + inverse).astype(np.uint32)
        vertices = np.concatenate((vertices, vertices[endpoints].mean(axis=1)))
        normals = np.concatenate((normals, normals[endpoints].mean(axis=1)))
        colors = np.concatenate((colors, colors[endpoints].mean(axis=1)))
        owners = np.concatenate((owners, owners[endpoints[:, 0]]))
        faces = np.vstack((np.c_[a, mid, c], np.c_[mid, b, c]))
    else:
        raise ValueError("Material sampling did not converge near the clipping plane")
    if len(faces):
        finished.append(faces)
    faces = np.concatenate(finished) if finished else np.empty((0, 3), np.uint32)
    used, indices = np.unique(faces, return_inverse=True)
    return Mesh(
        vertices[used],
        normals[used],
        colors[used],
        indices.reshape(-1, 3),
        owners[used],
        mesh.opacity,
    )


def contour_colors(colors, mesh, matrices, edges, width, color, background, fog_range):
    """Retain contour ink where native transparent bodies cover GPU contours."""
    from scipy.spatial import cKDTree

    from .export import visible_edges

    if not len(colors):
        return colors
    orthoscopic = abs(matrices[1][3, 3]) > 0.5
    edges = visible_edges(edges, matrices[0], orthoscopic, False)
    if not len(edges):
        return colors
    points, depth, eye = project(edges[:, :6].reshape(-1, 3), matrices)
    points, depth, eye = (a.reshape(-1, 2, a.shape[-1]) for a in (points, depth, eye))
    keep = (eye[:, :, 2] < 0).all(axis=1) & (np.abs(depth[:, :, 2]) <= 1).all(axis=1)
    points, depth, eye = points[keep], depth[keep], eye[keep]
    if not len(points):
        return colors
    lengths = np.linalg.norm(points[:, 1] - points[:, 0], axis=1)
    counts = np.maximum(np.ceil(lengths * 3).astype(int), 1)
    if counts.sum() > 10000000:
        raise ValueError("Projected contours exceed the material sampling budget")
    owners = np.repeat(np.arange(len(points)), counts)
    fractions = np.concatenate([(np.arange(n) + 0.5) / n for n in counts])
    tree = cKDTree(
        points[owners, 0] + fractions[:, None] * (points[owners, 1] - points[owners, 0])
    )
    pixels, ndc, positions = project(mesh.vertices, matrices)
    scale = 2 / (matrices[1][1, 1] * matrices[2][3])
    for start in range(0, len(pixels), 65536):
        stop = start + 65536
        xy = pixels[start:stop]
        _, closest = tree.query(xy, k=min(4, len(owners)))
        closest = np.asarray(closest).reshape(len(xy), -1)
        ids = owners[closest]
        a, b = points[ids, 0], points[ids, 1]
        delta = b - a
        t = np.clip(
            np.sum((xy[:, None] - a) * delta, axis=2)
            / np.maximum(np.sum(delta**2, axis=2), 1e-12),
            0,
            1,
        )
        distance = np.linalg.norm(xy[:, None] - a - t[:, :, None] * delta, axis=2)
        z = depth[ids, 0, 2] + t * (depth[ids, 1, 2] - depth[ids, 0, 2])
        eye_z = eye[ids, 0, 2] + t * (eye[ids, 1, 2] - eye[ids, 0, 2])
        pixel_size = scale * (1 if orthoscopic else -eye_z)
        radius = np.maximum(width / pixel_size, 1) / 2
        slope = (
            abs(matrices[1][2, 2])
            if orthoscopic
            else abs(matrices[1][2, 3]) / np.maximum(eye_z**2, 1e-8)
        )
        visible = (distance <= radius) & (
            z <= ndc[start:stop, None, 2] + np.maximum(width, pixel_size) * slope
        )
        ink = visible.any(axis=1)
        factor = np.clip(
            (fog_range[1] + positions[start:stop, 2])
            / max(fog_range[1] - fog_range[0], 1e-8),
            0,
            1,
        )
        colors[start:stop][ink] = (
            np.asarray(background)
            + (np.asarray(color) - background) * factor[ink, None]
        )
    return colors


def material_samples(
    mesh, material, matrices, background, budget, fog_range=(0, 1e10), contours=None
):
    mesh = clip_mesh(mesh, matrices)
    pixels, ndc, _ = project(mesh.vertices, matrices)
    original_points = np.c_[pixels * 3, ndc[:, 2]]
    original_faces = mesh.faces
    width, height = map(int, matrices[2][2:] * 3)
    if width * height * 12 > budget:
        raise ValueError("Visibility samples exceed cache_mb; reduce the image size")
    image, used = visibility(original_points, original_faces, width, height)
    mesh = Mesh(
        mesh.vertices,
        mesh.normals,
        mesh.colors,
        mesh.faces[used],
        mesh.owners,
        mesh.opacity,
    )
    mesh = subdivide(mesh, matrices, budget)
    pixels, ndc, _ = project(mesh.vertices, matrices)
    keep = visible_samples(
        np.c_[pixels * 3, ndc[:, 2]], mesh.faces, original_points, original_faces, image
    )
    mesh = Mesh(
        mesh.vertices,
        mesh.normals,
        mesh.colors,
        mesh.faces[keep],
        mesh.owners,
        mesh.opacity,
    )
    pixels, _, eye = project(mesh.vertices, matrices)
    rotation = matrices[0][:3, :3]
    direction = None if abs(matrices[1][3, 3]) > 0.5 else unit(-eye)
    fog = np.clip(
        (fog_range[1] + eye[:, 2]) / max(fog_range[1] - fog_range[0], 1e-8), 0, 1
    )
    if material == "richardson":
        tone = drawing_tone(mesh.normals @ rotation.T, direction, fog)
        colors = pencil(pixels, tone, mesh.colors, False)
    else:
        colors = bake(mesh, material, rotation, background, direction, fog)
    if contours is not None:
        colors = contour_colors(
            colors, mesh, matrices, *contours, background, fog_range
        )
    if mesh.opacity < 0.999999:
        mesh.opacity, colors = native_blend(colors, mesh.opacity, background)
    return mesh, colors


def contour_bands(
    mesh, edges, matrices, width, color, background, fog_range, budget, outer_only=False
):
    """Project visible contour bands with the same depth rule as the GPU pass."""
    from .export import visible_edges

    orthoscopic = abs(matrices[1][3, 3]) > 0.5
    edges = visible_edges(edges, matrices[0], orthoscopic, False)
    if not len(edges):
        return Mesh([], [], [], [], [])
    transform = matrices[1] @ matrices[0]
    endpoints = edges[:, :6].reshape(-1, 2, 3).astype(float)
    clip = (
        np.concatenate((endpoints, np.ones((*endpoints.shape[:2], 1))), axis=2)
        @ transform.T
    )
    lower, upper = np.zeros(len(edges)), np.ones(len(edges))
    for axis in range(3):
        for sign in (-1, 1):
            distance = clip[:, :, 3] + sign * clip[:, :, axis]
            delta = distance[:, 1] - distance[:, 0]
            crossing = np.divide(
                -distance[:, 0], delta, out=np.zeros_like(delta), where=delta != 0
            )
            lower = np.where(
                (distance[:, 0] < 0) & (delta > 0), np.maximum(lower, crossing), lower
            )
            upper = np.where(
                (distance[:, 1] < 0) & (delta < 0), np.minimum(upper, crossing), upper
            )
            upper[(distance < 0).all(axis=1)] = -1
    keep = lower <= upper
    endpoints, lower, upper = endpoints[keep], lower[keep], upper[keep]
    delta = endpoints[:, 1] - endpoints[:, 0]
    endpoints = (
        endpoints[:, :1] + np.stack((lower, upper), axis=1)[:, :, None] * delta[:, None]
    )
    xy, ndc, _ = project(endpoints.reshape(-1, 3), matrices)
    xy, ndc = xy.reshape(-1, 2, 2), ndc.reshape(-1, 2, 3)
    count = np.maximum(
        1, np.ceil(np.linalg.norm(xy[:, 1] - xy[:, 0], axis=1) * 3).astype(int)
    )
    if count.sum() * 600 > budget:
        raise ValueError("Ray contour bands exceed cache_mb")
    ids = np.repeat(np.arange(len(xy)), count)
    t = np.concatenate([np.arange(n) / n for n in count]) if len(count) else np.empty(0)
    t = np.stack((t, t + 1 / count[ids]), axis=1)
    positions = xy[ids, :1] + t[:, :, None] * (xy[ids, 1:] - xy[ids, :1])
    z = ndc[ids, 0, 2, None] + t * (ndc[ids, 1, 2] - ndc[ids, 0, 2])[:, None]
    clipped = clip_mesh(mesh, matrices)
    pixels, depth, _ = project(clipped.vertices, matrices)
    w, h = map(int, matrices[2][2:] * 3)
    if w * h * 12 > budget:
        raise ValueError("Ray contour visibility exceeds cache_mb")
    _, _, depth = visibility(np.c_[pixels * 3, depth[:, 2]], clipped.faces, w, h, True)

    def lookup(points, dx=0, dy=0):
        x = np.clip(np.floor(points[..., 0] * 3).astype(int) + dx, 0, w - 1)
        y = np.clip(np.floor(points[..., 1] * 3).astype(int) + dy, 0, h - 1)
        return depth[y, x]

    middle = positions.mean(axis=1)
    farthest = np.maximum.reduce(
        [lookup(middle, dx, dy) for dy in (-1, 0, 1) for dx in (-1, 0, 1)]
    )
    keep = z.mean(axis=1) <= farthest + 4e-6
    if outer_only:
        eye_z = -(fog_range[0] + 0.2 * (fog_range[1] - fog_range[0]))
        limit = (matrices[1][2, 2] * eye_z + matrices[1][2, 3]) / (
            matrices[1][3, 2] * eye_z + matrices[1][3, 3]
        )
        keep &= farthest >= limit
    positions, z = positions[keep], z[keep]
    eye_z = np.full_like(z, -1)
    if not orthoscopic:
        eye_z = -matrices[1][2, 3] / (z + matrices[1][2, 2])
    else:
        eye_z = (z - matrices[1][2, 3]) / matrices[1][2, 2]
    scale = matrices[1][1, 1] * matrices[2][3] / (2 * (1 if orthoscopic else -eye_z))
    radius = np.broadcast_to(np.maximum(width * scale, 1) / 2, z.shape)
    tangent = unit(positions[:, 1] - positions[:, 0])
    side = np.c_[-tangent[:, 1], tangent[:, 0]]
    offsets = side[:, None] * radius[:, :, None]
    corners = np.stack(
        (
            positions[:, 0] - offsets[:, 0],
            positions[:, 0] + offsets[:, 0],
            positions[:, 1] - offsets[:, 1],
            positions[:, 1] + offsets[:, 1],
        ),
        axis=1,
    )
    cz = np.minimum(z[:, [0, 0, 1, 1]], lookup(corners)) - 4e-7
    clip = np.concatenate(
        (
            np.stack(
                (
                    2 * corners[:, :, 0] / matrices[2][2] - 1,
                    1 - 2 * corners[:, :, 1] / matrices[2][3],
                    cz,
                ),
                axis=2,
            ),
            np.ones((*cz.shape, 1)),
        ),
        axis=2,
    )
    vertices = clip @ np.linalg.inv(transform).T
    vertices = (vertices[:, :, :3] / vertices[:, :, 3:]).reshape(-1, 3)
    fog = np.clip((fog_range[1] + eye_z) / max(fog_range[1] - fog_range[0], 1e-8), 0, 1)
    colors = np.asarray(background) + (np.asarray(color) - background) * fog[:, :, None]
    colors = colors[:, [0, 0, 1, 1]].reshape(-1, 3)
    a = np.arange(len(corners)) * 4
    faces = np.vstack((np.c_[a, a + 1, a + 2], np.c_[a + 1, a + 3, a + 2]))
    return Mesh(
        vertices,
        np.tile(matrices[0][2, :3], (len(vertices), 1)),
        colors,
        faces,
        np.zeros(len(vertices), int),
    )
