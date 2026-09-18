"""Managed views with reversible native state and session restoration."""

import re
from dataclasses import dataclass, field
from time import perf_counter

import numpy as np

from . import export, geometry, source
from .gpu import Drawing, Piece, Pool
from .presets import COLORS, QUALITIES, resolve


@dataclass
class Entry:
    name: str
    options: dict
    saved: dict
    drawings: dict
    object_settings: dict
    generated: list = field(default_factory=list)
    seconds: float = 0.0
    sampling_seconds: float = 0.0

    @property
    def nbytes(self):
        return sum(
            p.mesh.nbytes for ds in self.drawings.values() for d in ds for p in d.pieces
        )

    @property
    def cgo_nbytes(self):
        return sum(export.native_cgo_bytes(ds) for ds in self.drawings.values())


class Manager:
    def __init__(self, cmd):
        self.cmd = cmd
        self.entries = {}
        self.pool = Pool()
        self.widget = None
        self.events = None
        self.busy = False
        self.selected_keys = None

    def attach(self):
        if self.events is None:
            from .picking import attach

            self.widget, self.events = attach(self)

    def source_settings(self, obj, enabled):
        from pymol.setting import _get_index

        explicit = {row[0]: row[2] for row in self.cmd.get_object_settings(obj) or []}
        return {
            "state": explicit.get(_get_index("state")),
            "all_states": explicit.get(_get_index("all_states")),
            "enabled": obj in enabled,
        }

    def release_gpu(self):
        if self.pool.context is not None and self.widget is not None:
            self.cmd._call_with_opengl_context(self.pool.clear)

    def load_entry(self, entry):
        cmd = self.cmd
        for obj_index, (source_name, drawings) in enumerate(entry.drawings.items(), 1):
            name = f"{entry.name}_shape_{obj_index}"
            body = f"{entry.name}_alpha_{obj_index}"
            proxy = f"{entry.name}_ray_{obj_index}"
            has_alpha = any(
                p.mesh.opacity < 0.999999 for d in drawings for p in d.pieces
            )
            has_opaque = any(
                p.mesh.opacity >= 0.999999 for d in drawings for p in d.pieces
            )
            targets = [name]
            if has_alpha:
                targets.append(body)
            if has_opaque:
                targets.append(proxy)
            entry.generated.extend(targets)
            for state, drawing in enumerate(drawings, 1):
                drawing.name, drawing.state = name, state
                cmd.load_callback(drawing, name, state, 1, 0, 1, 0)
                if has_alpha:
                    values = []
                    for piece in drawing.pieces:
                        if piece.mesh.opacity < 0.999999:
                            values.extend(
                                export.cgo_mesh(piece, drawing.profile.material)
                            )
                    cmd.load_cgo(values, body, state=state, zoom=0)
                if has_opaque:
                    cmd.load_cgo(export.ray_proxy(drawing), proxy, state=state, zoom=0)
            for target in targets:
                cmd.group(entry.name, target)
                for setting, value in entry.object_settings[source_name].items():
                    if setting != "enabled" and value is not None:
                        cmd.set(setting, value, target)
                if not entry.object_settings[source_name]["enabled"]:
                    cmd.disable(target)
            if has_alpha:
                cmd.set("cgo_lighting", 0, body)
                cmd.set("cgo_transparency", 0, body)

    def unload_entry(self, entry):
        for name in entry.generated:
            self.cmd.delete(name)
        self.cmd.delete(entry.name)
        entry.generated.clear()

    def apply(
        self,
        style,
        selection,
        representation,
        color,
        quality,
        name,
        edge,
        edge_width,
        edge_color,
        transparency,
        cache_mb=2048,
    ):
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
            raise ValueError(
                "name must start with a letter and contain only letters, digits, and underscores"
            )
        if self.cmd.get_legal_name(name) != name:
            raise ValueError(
                "name is a reserved PyMOL word; choose another managed group name"
            )
        if color not in COLORS or quality not in QUALITIES:
            raise ValueError(
                f"color must be one of {COLORS}; quality must be one of {tuple(QUALITIES)}"
            )
        profile = resolve(style, representation, edge, edge_width)
        cache_mb = float(cache_mb)
        if not np.isfinite(cache_mb) or cache_mb <= 0:
            raise ValueError("cache_mb must be finite and positive")
        budget = int(cache_mb * 1024**2)
        if transparency is not None:
            transparency = float(transparency)
            if not np.isfinite(transparency) or not 0 <= transparency <= 1:
                raise ValueError("transparency must be between 0 and 1")
        edge_rgb = tuple(self.cmd.get_color_tuple(edge_color))
        options = {
            "style": style,
            "selection": selection,
            "representation": representation,
            "color": color,
            "quality": quality,
            "name": name,
            "edge": edge,
            "edge_width": edge_width,
            "edge_color": edge_color,
            "transparency": transparency,
            "cache_mb": cache_mb,
        }
        old = self.entries.get(name)
        if old is None and name in self.cmd.get_names("all"):
            raise ValueError(f"An object or selection named {name!r} already exists")
        self.busy = True
        start = perf_counter()
        candidate = None
        try:
            with export.settings(self.cmd, {"suspend_updates": 1}):
                if old:
                    source.restore_reps(self.cmd, old.saved)
                states, saved = source.read(self.cmd, selection, budget)
                for other_name, entry in self.entries.items():
                    if other_name != name and saved.keys() & entry.saved.keys():
                        raise ValueError(
                            f"Selection overlaps managed view {other_name!r}; reset it or use a disjoint selection"
                        )
                drawings = {}
                cache_size = 0
                for obj, frames in states.items():
                    drawings[obj] = []
                    for frame in frames:
                        pieces = []
                        for rep, subset, opacity in source.layers(
                            frame, profile.representation, transparency
                        ):
                            mesh = geometry.build(
                                subset.atoms,
                                subset.coords,
                                subset.bonds,
                                subset.model,
                                profile,
                                rep,
                                quality,
                                color,
                            )
                            mesh.opacity = opacity
                            if len(mesh.faces):
                                pieces.append(Piece(mesh, subset.atoms))
                                cache_size += mesh.nbytes
                                if cache_size > budget:
                                    raise ValueError(
                                        "Prepared geometry exceeds cache_mb; lower quality or increase cache_mb"
                                    )
                        drawings[obj].append(
                            Drawing(
                                pieces,
                                profile,
                                edge_rgb,
                                self.pool,
                                sample_budget=budget,
                                anchors=frame.coords,
                                keys=tuple((a.model, a.index) for a in frame.atoms),
                                background=tuple(
                                    self.cmd.get_color_tuple(self.cmd.get("bg_rgb"))
                                ),
                            )
                        )
                if (
                    not any(d.pieces for ds in drawings.values() for d in ds)
                    and transparency != 1
                ):
                    raise ValueError("The selected atoms produced no drawable geometry")
                enabled = set(self.cmd.get_names("objects", enabled_only=1))
                object_settings = {
                    obj: self.source_settings(obj, enabled) for obj in drawings
                }
                candidate = Entry(name, options, saved, drawings, object_settings)
                if candidate.cgo_nbytes > budget:
                    raise ValueError(
                        "Prepared native CGO exceeds cache_mb; lower quality or increase cache_mb"
                    )
                existing = set(self.cmd.get_names("all"))
                owned = set(old.generated) if old else set()
                for index, frames in enumerate(drawings.values(), 1):
                    generated = [f"{name}_shape_{index}"]
                    if any(p.mesh.opacity < 0.999999 for d in frames for p in d.pieces):
                        generated.append(f"{name}_alpha_{index}")
                    if any(
                        p.mesh.opacity >= 0.999999 for d in frames for p in d.pieces
                    ):
                        generated.append(f"{name}_ray_{index}")
                    if set(generated) & (existing - owned):
                        raise ValueError(
                            "A generated object name already exists; choose another name"
                        )
                if old:
                    self.unload_entry(old)
                self.load_entry(candidate)
                source.hide_reps(self.cmd, saved)
                self.entries[name] = candidate
                candidate.seconds = perf_counter() - start
        except Exception:
            if candidate:
                if candidate.generated:
                    self.unload_entry(candidate)
                source.restore_reps(self.cmd, candidate.saved)
            if old:
                if not old.generated:
                    self.load_entry(old)
                source.hide_reps(self.cmd, old.saved)
                self.entries[name] = old
            raise
        finally:
            self.busy = False
        self.release_gpu()
        self.attach()
        self.selected_keys = None
        self.update_selection()
        self.cmd.refresh()
        return candidate

    def reset(self, name="all", restore=True):
        names = list(self.entries) if name == "all" else [name]
        self.busy = True
        try:
            with export.settings(self.cmd, {"suspend_updates": 1}):
                for key in names:
                    entry = self.entries.pop(key, None)
                    if entry is not None:
                        self.unload_entry(entry)
                        if restore:
                            source.restore_reps(self.cmd, entry.saved)
            self.release_gpu()
            if not self.entries and self.events is not None:
                self.events.close()
                self.events = self.widget = None
        finally:
            self.busy = False
        self.cmd.refresh()

    def refresh(self, name):
        names = list(self.entries) if name == "all" else [name]
        for key in names:
            if key not in self.entries:
                raise ValueError(f"No managed view named {key!r}")
            self.apply(**self.entries[key].options)

    def active_drawings(self):
        enabled = set(self.cmd.get_names("objects", enabled_only=1))
        for entry in self.entries.values():
            if entry.name not in enabled:
                continue
            for source_name, drawings in entry.drawings.items():
                name = drawings[0].name
                if name not in enabled or source_name not in enabled:
                    continue
                if self.cmd.get_setting_int("all_states", name):
                    yield from drawings
                else:
                    state = self.cmd.get_object_state(name)
                    if len(drawings) == 1:
                        state = 1
                    if 1 <= state <= len(drawings):
                        yield drawings[state - 1]

    def maintenance(self):
        if self.busy:
            return
        names = set(self.cmd.get_names("objects"))
        enabled = set(self.cmd.get_names("objects", enabled_only=1))
        for key, entry in list(self.entries.items()):
            errors = [d.error for ds in entry.drawings.values() for d in ds if d.error]
            if errors:
                self.reset(key)
                print(
                    f" cuemol_style: restored native view {key!r} after a drawing failure: {errors[0]}"
                )
                continue
            if (
                key not in names
                or any(obj not in names for obj in entry.drawings)
                or any(n not in names for n in entry.generated)
            ):
                self.reset(key)
                continue
            for i, obj in enumerate(entry.drawings, 1):
                current = self.source_settings(obj, enabled)
                previous = entry.object_settings[obj]
                if current != previous:
                    for target in (
                        f"{key}_shape_{i}",
                        f"{key}_alpha_{i}",
                        f"{key}_ray_{i}",
                    ):
                        if target not in names:
                            continue
                        for setting in ("state", "all_states"):
                            if current[setting] != previous[setting]:
                                if current[setting] is None:
                                    self.cmd.unset(setting, target)
                                else:
                                    self.cmd.set(setting, current[setting], target)
                        if current["enabled"] != previous["enabled"]:
                            (
                                self.cmd.enable
                                if current["enabled"]
                                else self.cmd.disable
                            )(target)
                    entry.object_settings[obj] = current
        try:
            self.prepare_view()
        except Exception as exc:  # noqa: BLE001 - GUI callback boundary.
            for drawing in self.active_drawings():
                drawing.error = str(exc)
        self.update_selection()

    def prepare_view(self, width=0, height=0, force=False, native=True):
        if self.busy or self.widget is None:
            return
        from .raster import layers
        from .sampling import camera, clipping

        matrices = camera(self.cmd, width, height)
        background = tuple(self.cmd.get_color_tuple(self.cmd.get("bg_rgb")))
        view = self.cmd.get_view()
        fog = (-view[11], -view[11] + (view[16] - view[15]) / 2)
        native_fog = None
        density = self.cmd.get_setting_float("fog")
        if self.cmd.get_setting_int("depth_cue") and density != 0:
            near, far = clipping(view)
            start = near + (far - near) * self.cmd.get_setting_float("fog_start")
            native_fog = (
                start,
                start + (far - start) / density if density > 1e-8 else far,
            )
        key = tuple(array.tobytes() for array in matrices) + (background, native_fog)
        active = {id(d) for d in self.active_drawings()}
        self.busy = True
        try:
            for entry in self.entries.values():
                started = perf_counter()
                budget = int(entry.options["cache_mb"] * 1024**2)
                for index, drawings in enumerate(entry.drawings.values(), 1):
                    for drawing in drawings:
                        drawing.background = background
                        drawing.fog = fog
                        drawing.sample_matrices = matrices
                        drawing.sample_budget = budget
                        if not native:
                            continue
                        if id(drawing) not in active and drawing.sampled_bytes:
                            values = []
                            for piece in drawing.pieces:
                                if piece.mesh.opacity < 0.999999:
                                    values.extend(
                                        export.cgo_mesh(piece, drawing.profile.material)
                                    )
                            self.cmd.load_cgo(
                                values,
                                f"{entry.name}_alpha_{index}",
                                state=drawing.state,
                                zoom=0,
                            )
                            self.cmd.set(
                                "cgo_lighting", 0, f"{entry.name}_alpha_{index}"
                            )
                            drawing.sampled_bytes = 0
                            drawing.sample_key = None
                        if id(drawing) not in active or (
                            not force and drawing.sample_key == key
                        ):
                            continue
                        values = []
                        for opacity, meshes in layers(drawing).items():
                            if opacity < 0.999999:
                                values.extend(
                                    export.sampled_cgo(
                                        drawing,
                                        meshes,
                                        matrices,
                                        budget,
                                        opacity,
                                        native_fog,
                                    )
                                )
                        if values:
                            size = len(values) * 4
                            base_alpha = sum(
                                4 * (5 + 36 * len(p.mesh.faces))
                                for p in drawing.pieces
                                if p.mesh.opacity < 0.999999
                            )
                            if (
                                entry.cgo_nbytes
                                - max(drawing.sampled_bytes, base_alpha)
                                + size
                                > budget
                            ):
                                raise ValueError(
                                    "Camera-dependent CGO exceeds cache_mb"
                                )
                            self.cmd.load_cgo(
                                values,
                                f"{entry.name}_alpha_{index}",
                                state=drawing.state,
                                zoom=0,
                            )
                            self.cmd.set(
                                "cgo_lighting", 0, f"{entry.name}_alpha_{index}"
                            )
                            drawing.sampled_bytes = size
                            entry.sampling_seconds = perf_counter() - started
                        drawing.sample_key = key
        finally:
            self.busy = False

    def update_selection(self):
        names = [
            name
            for name in self.cmd.get_names("selections", enabled_only=1)
            if not name.startswith("_")
        ]
        keys = (
            frozenset(self.cmd.index(" or ".join("%" + name for name in names)))
            if names
            else frozenset()
        )
        if keys == self.selected_keys:
            return
        self.selected_keys = keys
        for entry in self.entries.values():
            for drawings in entry.drawings.values():
                for drawing in drawings:
                    indices = [i for i, key in enumerate(drawing.keys) if key in keys]
                    drawing.selected = np.ascontiguousarray(drawing.anchors[indices])
        self.cmd.refresh()


def manager_for(cmd):
    owner = cmd._pymol
    manager = vars(owner).get("_cuemol_style_manager")
    if manager is None:
        manager = owner._cuemol_style_manager = Manager(cmd)
    for attr, hook in (
        ("_session_save_tasks", save_session),
        ("_session_restore_tasks", restore_session),
    ):
        tasks = getattr(owner, attr)
        if hook not in tasks:
            tasks.append(hook)
    return manager


def save_session(session, _self):
    manager = manager_for(_self)
    session["cuemol_style"] = {
        "version": 2,
        "entries": [
            {"options": e.options, "saved": e.saved, "generated": e.generated}
            for e in manager.entries.values()
        ],
    }
    return 1


def restore_session(session, _self):
    manager = manager_for(_self)
    if manager.events is not None:
        manager.events.close()
    manager.release_gpu()
    manager.entries.clear()
    manager.widget = manager.events = None
    data = session.get("cuemol_style", {})
    for entry in data.get("entries", []):
        for name in entry["generated"]:
            _self.delete(name)
        _self.delete(entry["options"]["name"])
        source.restore_reps(_self, entry["saved"])
    for entry in data.get("entries", []):
        manager.apply(**entry["options"])
    return 1
