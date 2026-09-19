"""Compact analytic atoms for large scenes, native ray, and source picking."""

from dataclasses import dataclass

import numpy as np


def transparency(cmd, owner, needed):
    """Share and restore PyMOL's shader-capable transparency mode across views."""
    session = cmd._pymol.session
    key = "_pymol_style_native_transparency"
    state = getattr(session, key, None)
    if needed:
        if state is None:
            state = {"previous": cmd.get_setting_int("transparency_mode"), "owners": []}
            setattr(session, key, state)
        if owner not in state["owners"]:
            state["owners"].append(owner)
        if cmd.get_setting_int("transparency_mode") != 3:
            cmd.set("transparency_mode", 3)
    elif state is not None:
        state["owners"] = [value for value in state["owners"] if value != owner]
        if not state["owners"]:
            if cmd.get_setting_int("transparency_mode") == 3:
                cmd.set("transparency_mode", state["previous"])
            delattr(session, key)


def groups(drawings, prefix):
    """Keep each native CGO at one object-level opacity for shader rendering."""
    opacities = sorted(
        {
            part.opacity
            for drawing in drawings
            for part in drawing.native
            if part.opacity > 0
        }
    )
    return {f"{prefix}_{index}": alpha for index, alpha in enumerate(opacities)}


@dataclass
class NativeAtoms:
    atoms: tuple
    spheres: np.ndarray
    sphere_colors: np.ndarray
    sphere_owners: np.ndarray
    cylinders: np.ndarray
    cylinder_colors: np.ndarray
    cylinder_owners: np.ndarray
    opacity: float = 1.0

    def __post_init__(self):
        if not np.isfinite(self.spheres).all() or not np.isfinite(self.cylinders).all():
            raise ValueError("Native atom coordinates and radii must be finite")
        if np.any(self.spheres[:, 3] < 0) or np.any(self.cylinders[:, 6] < 0):
            raise ValueError("Native atom radii must be non-negative")

    @property
    def nbytes(self):
        return sum(
            a.nbytes
            for a in (
                self.spheres,
                self.sphere_colors,
                self.sphere_owners,
                self.cylinders,
                self.cylinder_colors,
                self.cylinder_owners,
            )
        )

    @property
    def cgo_nbytes(self):
        return 4 * (2 + 9 * len(self.spheres) + 14 * len(self.cylinders))

    @property
    def extent(self):
        low, high = [], []
        if len(self.spheres):
            low.append((self.spheres[:, :3] - self.spheres[:, 3:4]).min(axis=0))
            high.append((self.spheres[:, :3] + self.spheres[:, 3:4]).max(axis=0))
        if len(self.cylinders):
            a, b, radius = (
                self.cylinders[:, :3],
                self.cylinders[:, 3:6],
                self.cylinders[:, 6:7],
            )
            low.append((np.minimum(a, b) - radius).min(axis=0))
            high.append((np.maximum(a, b) + radius).max(axis=0))
        return [np.min(low, axis=0), np.max(high, axis=0)] if low else [[0] * 3] * 2

    def cgo(self, opacity=None):
        from pymol.cgo import ALPHA, COLOR, CYLINDER, SPHERE

        spheres = np.empty((len(self.spheres), 9), np.float32)
        spheres[:, 0], spheres[:, 4] = COLOR, SPHERE
        spheres[:, 1:4], spheres[:, 5:9] = self.sphere_colors, self.spheres
        cylinders = np.empty((len(self.cylinders), 14), np.float32)
        cylinders[:, 0] = CYLINDER
        cylinders[:, 1:8] = self.cylinders
        cylinders[:, 8:11] = cylinders[:, 11:14] = self.cylinder_colors
        return np.concatenate(
            (
                [ALPHA, self.opacity if opacity is None else opacity],
                spheres.ravel(),
                cylinders.ravel(),
            )
        ).tolist()

    def hit(self, origin, direction):
        """Intersect exact spheres and capped cylinders without tessellation."""
        best = None

        def retain(distances, owners):
            nonlocal best
            if not len(distances):
                return
            distances = np.where(distances > 1e-8, distances, np.inf)
            index = int(np.argmin(distances))
            distance = float(distances[index])
            if np.isfinite(distance) and (best is None or distance < best[0]):
                best = distance, int(owners[index])

        if self.opacity <= 0:
            return None
        if len(self.spheres):
            delta = np.asarray(origin) - self.spheres[:, :3]
            b = delta @ direction
            discriminant = (
                b * b - np.sum(delta * delta, axis=1) + self.spheres[:, 3] ** 2
            )
            root = np.sqrt(np.maximum(discriminant, 0))
            distance = np.where(-b - root > 1e-8, -b - root, -b + root)
            retain(np.where(discriminant >= 0, distance, np.inf), self.sphere_owners)
        if len(self.cylinders):
            start = self.cylinders[:, :3].astype(float)
            axis = self.cylinders[:, 3:6] - start
            length = np.linalg.norm(axis, axis=1)
            axis /= np.maximum(length[:, None], 1e-12)
            delta = np.asarray(origin) - start
            axial = axis @ direction
            position = np.sum(delta * axis, axis=1)
            a = 1 - axial * axial
            b = delta @ direction - position * axial
            c = (
                np.sum(delta * delta, axis=1)
                - position * position
                - self.cylinders[:, 6] ** 2
            )
            discriminant = b * b - a * c
            root = np.sqrt(np.maximum(discriminant, 0))
            for sign in (-1, 1):
                distance = (-b + sign * root) / np.maximum(a, 1e-12)
                height = position + distance * axial
                valid = (
                    (a > 1e-12)
                    & (discriminant >= 0)
                    & (height >= 0)
                    & (height <= length)
                )
                retain(np.where(valid, distance, np.inf), self.cylinder_owners)
            for height in (np.zeros_like(length), length):
                distance = np.divide(
                    height - position,
                    axial,
                    out=np.zeros_like(length),
                    where=abs(axial) > 1e-12,
                )
                point = delta + distance[:, None] * direction - height[:, None] * axis
                valid = (abs(axial) > 1e-12) & (
                    np.sum(point * point, axis=1) <= self.cylinders[:, 6] ** 2
                )
                retain(np.where(valid, distance, np.inf), self.cylinder_owners)
        return best


class Builder:
    """Collect primitive parameters without allocating individual meshes."""

    def __init__(self, atoms):
        self.atoms = tuple(atoms)
        self.spheres, self.sphere_colors, self.sphere_owners = [], [], []
        self.cylinders, self.cylinder_colors, self.cylinder_owners = [], [], []

    def sphere(self, center, radius, color, owner, detail):
        self.spheres.append((*center, radius))
        self.sphere_colors.append(color)
        self.sphere_owners.append(owner)

    def cylinder(self, a, b, radius, color, owner=0, detail=12):
        if np.linalg.norm(np.asarray(b) - a) < 1e-9:
            return self.sphere(a, radius, color, owner, detail)
        self.cylinders.append((*a, *b, radius))
        self.cylinder_colors.append(color)
        self.cylinder_owners.append(owner)

    def dashed(self, a, b, radius, color, owner=0, detail=8, count=7):
        for index in range(count):
            self.cylinder(
                a + (b - a) * index / count,
                a + (b - a) * (index + 0.5) / count,
                radius,
                color,
                owner,
                detail,
            )

    def finish(self):
        return NativeAtoms(
            self.atoms,
            np.asarray(self.spheres, np.float32).reshape(-1, 4),
            np.asarray(self.sphere_colors, np.float32).reshape(-1, 3),
            np.asarray(self.sphere_owners, np.int32),
            np.asarray(self.cylinders, np.float32).reshape(-1, 7),
            np.asarray(self.cylinder_colors, np.float32).reshape(-1, 3),
            np.asarray(self.cylinder_owners, np.int32),
        )
