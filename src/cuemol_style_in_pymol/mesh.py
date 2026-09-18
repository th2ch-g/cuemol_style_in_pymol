"""Indexed geometry shared by the GPU renderer, picker, and CGO exporter."""

from dataclasses import dataclass
from functools import cached_property

import numpy as np


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-12)


@dataclass
class Mesh:
    vertices: np.ndarray
    normals: np.ndarray
    colors: np.ndarray
    faces: np.ndarray
    owners: np.ndarray
    opacity: float = 1.0

    def __post_init__(self):
        self.vertices = np.ascontiguousarray(self.vertices, dtype=np.float32).reshape(
            -1, 3
        )
        self.normals = np.ascontiguousarray(self.normals, dtype=np.float32).reshape(
            -1, 3
        )
        self.colors = np.ascontiguousarray(self.colors, dtype=np.float32).reshape(-1, 3)
        self.faces = np.ascontiguousarray(self.faces, dtype=np.uint32).reshape(-1, 3)
        self.owners = np.asarray(self.owners, dtype=np.int32)
        if not np.isfinite(self.vertices).all() or not np.isfinite(self.normals).all():
            raise ValueError("Non-finite mesh coordinates or normals")

    @property
    def nbytes(self):
        return sum(
            a.nbytes
            for a in (self.vertices, self.normals, self.colors, self.faces, self.owners)
        ) + sum(a.nbytes for a in self.__dict__.get("face_bounds", ()))

    @cached_property
    def face_bounds(self):
        """Conservative boxes for contiguous batches shared by culling and picking."""
        if not len(self.faces):
            return np.empty((0, 3), np.float32), np.empty((0, 3), np.float32)
        triangles = self.vertices[self.faces]
        starts = np.arange(0, len(triangles), 64)
        return (
            np.minimum.reduceat(triangles.min(axis=1), starts),
            np.maximum.reduceat(triangles.max(axis=1), starts),
        )

    def edge_data(self):
        """Weld section seams before finding boundary and crease edges."""
        _, welded = np.unique(np.round(self.vertices, 5), axis=0, return_inverse=True)
        faces = self.faces.astype(int)
        area = np.cross(
            self.vertices[faces[:, 1]] - self.vertices[faces[:, 0]],
            self.vertices[faces[:, 2]] - self.vertices[faces[:, 0]],
        )
        faces = faces[np.linalg.norm(area, axis=1) > 1e-8]
        normals = unit(
            np.cross(
                self.vertices[faces[:, 1]] - self.vertices[faces[:, 0]],
                self.vertices[faces[:, 2]] - self.vertices[faces[:, 0]],
            )
        )
        if not len(faces):
            return np.empty((0, 12), np.float32)
        edges = np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2)
        keys = np.sort(welded[edges], axis=1)
        order = np.lexsort((keys[:, 1], keys[:, 0]))
        keys, edges = keys[order], edges[order]
        fn = np.repeat(normals, 3, axis=0)[order]
        first = np.r_[0, np.flatnonzero(np.any(np.diff(keys, axis=0), axis=1)) + 1]
        last = np.r_[first[1:] - 1, len(keys) - 1]
        n0 = fn[first]
        n1 = np.where((last == first)[:, None], -n0, fn[last])
        keep = (np.sum(n0 * n1, axis=1) < 0.9999) & (keys[first, 0] != keys[first, 1])
        return np.ascontiguousarray(
            np.c_[
                self.vertices[edges[first, 0]], self.vertices[edges[first, 1]], n0, n1
            ][keep],
            dtype=np.float32,
        )


def merge(meshes, opacity=1.0):
    meshes = [m for m in meshes if len(m.faces)]
    if not meshes:
        return Mesh([], [], [], [], [], opacity)
    offsets = np.cumsum([0] + [len(m.vertices) for m in meshes[:-1]])
    return Mesh(
        np.concatenate([m.vertices for m in meshes]),
        np.concatenate([m.normals for m in meshes]),
        np.concatenate([m.colors for m in meshes]),
        np.concatenate([m.faces + o for m, o in zip(meshes, offsets)]),
        np.concatenate([m.owners for m in meshes]),
        opacity,
    )


def ray_hits(mesh, origin, direction):
    """Return the nearest positive Moller-Trumbore intersection and owner."""
    low, high = mesh.face_bounds
    parallel = np.abs(direction) < 1e-12
    inverse = np.divide(1.0, direction, out=np.ones(3), where=~parallel)
    a, b = (low - origin) * inverse, (high - origin) * inverse
    near, far = np.minimum(a, b), np.maximum(a, b)
    outside = (origin < low) | (origin > high)
    near[:, parallel] = np.where(outside[:, parallel], np.inf, -np.inf)
    far[:, parallel] = np.where(outside[:, parallel], -np.inf, np.inf)
    groups = np.flatnonzero(far.min(axis=1) >= np.maximum(0, near.max(axis=1)))
    if not len(groups):
        return None
    indices = (groups[:, None] * 64 + np.arange(64)).ravel()
    indices = indices[indices < len(mesh.faces)]
    p = mesh.vertices[mesh.faces[indices]]
    a, b = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    h = np.cross(direction, b)
    det = np.einsum("ij,ij->i", a, h)
    good = np.abs(det) > 1e-9
    inv = np.divide(1.0, det, out=np.zeros_like(det), where=good)
    s = origin - p[:, 0]
    u = inv * np.einsum("ij,ij->i", s, h)
    q = np.cross(s, a)
    v = inv * (q @ direction)
    t = inv * np.einsum("ij,ij->i", b, q)
    good &= (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
    if not good.any():
        return None
    index = np.argmin(np.where(good, t, np.inf))
    return float(t[index]), int(mesh.owners[mesh.faces[indices[index], 0]])
