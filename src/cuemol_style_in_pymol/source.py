"""Read-only molecular snapshots and reversible representation changes."""

from dataclasses import dataclass
from copy import deepcopy
from contextlib import contextmanager
from uuid import uuid4

import numpy as np


# PyMOL 3.x representation bits, independent of user-visible object types.
REPS = {"sticks": 1, "cpk": 2, "surface": 4, "cartoon": 32, "ribbon": 64}
REPLACED = sum(REPS.values()) | 16 | 128 | 2048


@dataclass(frozen=True, slots=True)
class Atom:
    model: str
    index: int
    name: str
    resi: str
    resn: str
    chain: str
    segi: str
    element: str
    ss: str
    kind: str
    color: tuple
    vdw: float
    alt: str = ""
    occupancy: float = 1.0


@dataclass
class State:
    atoms: tuple
    coords: np.ndarray
    bonds: tuple
    model: object
    representations: np.ndarray
    transparency: dict

    def subset(self, mask):
        indices = np.flatnonzero(mask)
        lookup = {int(old): new for new, old in enumerate(indices)}
        model = deepcopy(self.model)
        model.atom = [model.atom[i] for i in indices]
        model.bond = [b for b in model.bond if all(i in lookup for i in b.index)]
        for b in model.bond:
            b.index = [lookup[i] for i in b.index]
        return State(
            tuple(self.atoms[i] for i in indices),
            self.coords[indices],
            tuple(tuple(b.index) for b in model.bond),
            model,
            self.representations[indices],
            {key: value[indices] for key, value in self.transparency.items()},
        )


@contextmanager
def atom_selection(cmd, keys):
    """Build selections from integer indices, without interpolating expressions."""
    name = "_cuemol_select_" + uuid4().hex
    part = name + "_part"
    objects = {}
    for model, index in keys:
        objects.setdefault(model, []).append(index)
    try:
        cmd.select(name, "none", quiet=1)
        for model, indices in objects.items():
            cmd.select_list(part, model, indices, mode="index", quiet=1)
            cmd.select(name, part, merge=1, quiet=1)
        yield name
    finally:
        cmd.delete(name)
        cmd.delete(part)


def read(cmd, selection, budget_bytes=None):
    """Capture all states before entering any OpenGL callback."""
    objects = cmd.get_object_list(f"({selection})")
    if not objects:
        raise ValueError("The selection contains no molecular objects")
    rows = []
    cmd.iterate(
        selection,
        "out.append((model,index,reps,color,s.cartoon_transparency,"
        "s.stick_transparency,s.sphere_transparency,s.transparency))",
        space={"out": rows},
    )
    if not rows:
        raise ValueError("The selection contains no atoms")
    properties = {(r[0], r[1]): r[2:] for r in rows}
    protein = set(cmd.index(f"({selection}) and polymer.protein"))
    nucleic = set(cmd.index(f"({selection}) and polymer.nucleic"))
    colors = {r[3]: tuple(cmd.get_color_tuple(r[3])) for r in rows}
    result = {}
    atom_cache = {}
    estimated_bytes = 0
    for obj in objects:
        keys = [key for key in properties if key[0] == obj]
        count = cmd.count_states(obj)
        estimated_bytes += len(keys) * count * 400
        if budget_bytes is not None and estimated_bytes > budget_bytes:
            raise ValueError(
                "Molecular snapshots exceed cache_mb; increase cache_mb or select fewer atoms"
            )
        states = []
        with atom_selection(cmd, keys) as sele:
            for state in range(1, count + 1):
                model = cmd.get_model(sele, state=state)
                atoms, prop = [], []
                for a in model.atom:
                    key = (obj, a.index)
                    p = properties[key]
                    kind = (
                        "protein"
                        if key in protein
                        else "nucleic"
                        if key in nucleic
                        else "other"
                    )
                    if key not in atom_cache:
                        atom_cache[key] = Atom(
                            obj,
                            a.index,
                            a.name,
                            a.resi,
                            a.resn,
                            a.chain,
                            a.segi,
                            a.symbol,
                            a.ss,
                            kind,
                            colors[p[1]],
                            a.vdw,
                            a.alt,
                            a.q,
                        )
                    atoms.append(atom_cache[key])
                    prop.append(p)
                prop = np.asarray(prop, float).reshape(-1, 6)
                states.append(
                    State(
                        tuple(atoms),
                        np.asarray(model.get_coord_list(), dtype=np.float32).reshape(
                            -1, 3
                        ),
                        tuple(tuple(b.index) for b in model.bond),
                        model,
                        prop[:, 0].astype(int),
                        {
                            rep: prop[:, col]
                            for rep, col in (
                                ("ribbon", 2),
                                ("cartoon", 2),
                                ("tube", 2),
                                ("nucleic", 2),
                                ("sticks", 3),
                                ("ballstick", 3),
                                ("cpk", 4),
                                ("surface", 5),
                            )
                        },
                    )
                )
        result[obj] = states
    return result, {key: int(p[0]) for key, p in properties.items()}


def layers(state, representation, transparency):
    """Keep mixed native representations and their per-atom opacity."""
    if not state.atoms:
        return
    if representation == "auto":
        choices = [
            (rep, (state.representations & mask) != 0) for rep, mask in REPS.items()
        ]
        choices = [(rep, mask) for rep, mask in choices if mask.any()]
        covered = np.zeros(len(state.atoms), bool)
        for _, mask in choices:
            covered |= mask
        fallback = ~covered & ((state.representations & (16 | 128 | 2048)) != 0)
        if fallback.any():
            choices.append(("ribbon", fallback))
    else:
        choices = [(representation, np.ones(len(state.atoms), bool))]
    for rep, mask in choices:
        opacity = (
            1 - state.transparency[rep]
            if transparency is None
            else np.full(len(mask), 1 - transparency)
        )
        # Six decimal places avoid spurious groups from setting serialization.
        for alpha in np.unique(np.round(opacity[mask], 6)):
            if alpha <= 0:
                continue
            subset = state.subset(mask & np.isclose(opacity, alpha, atol=5e-7, rtol=0))
            yield rep, subset, float(alpha)


def restore_reps(cmd, saved):
    existing = set(cmd.get_object_list())
    keys = [key for key in saved if key[0] in existing]
    if keys:
        with atom_selection(cmd, keys) as sele:
            cmd.alter(
                sele,
                "reps=saved.get((model,index),reps)",
                space={"saved": saved},
                quiet=1,
            )
            cmd.rebuild(sele)


def hide_reps(cmd, saved):
    with atom_selection(cmd, saved) as sele:
        cmd.alter(sele, "reps=reps & ~mask", space={"mask": REPLACED}, quiet=1)
        cmd.rebuild(sele)
