"""CueMol-inspired molecular geometry and materials for PyMOL 3.1 Qt."""


def cuemol_style(
    style="ribbon",
    selection="all",
    representation=None,
    color="cuemol",
    quality="medium",
    name="cuemol",
    edge="auto",
    edge_width="",
    edge_color="black",
    transparency="keep",
    cache_mb=2048,
    filename="",
    width=0,
    height=0,
    quiet=0,
    _self=None,
):
    """
    DESCRIPTION

        Apply CueMol-inspired geometry and live GPU materials without changing
        PyMOL's source code or its standard commands. Requires PyMOL 3.1 Qt and
        compatibility OpenGL 2.1 / GLSL 1.20. CueMol itself is not required.

    USAGE

        cuemol_style style [, selection [, representation [, color [, quality [, name]]]]]
        cuemol_style list
        cuemol_style refresh [, name=cuemol]
        cuemol_style reset [, name=cuemol]
        cuemol_style png, filename=figure.png [, width=2400, height=1800]
        cuemol_style ray, filename=figure_ray.png [, width=2400, height=1800]

    ARGUMENTS

        style: ribbon (default), richardson, round_ribbon, fancy_ribbon, cartoon,
            round_cartoon, tube, nucleic, ballstick, cpk, surface, outline,
            silhouette, default, shadow, nolighting, matte, toon1, toon2,
            diff_metal, spec_metal, metallic_chrome, metallic_copper, stone35,
            wood31, wood14scl2. Use list to inspect profiles and active views.
        selection: molecular atoms to style {default: all}
        representation: auto, ribbon, cartoon, tube, nucleic, ballstick,
            sticks, cpk, surface {default: ribbon for material/outline styles;
            geometry presets select their own representation}
            auto inherits source layers for material/outline styles.
        color: cuemol, keep, chain, ss, rainbow, element {default: cuemol}
            cuemol follows CueMol GUI defaults: khaki helices, SteelBlue sheets,
            FloralWhite coils, yellow nucleic geometry, and DefaultCPKColoring
            for atomic representations (carbon inherits the molecular color).
            ss retains the optional WoodyHSCPaint palette.
        quality: low, medium, high {default: medium}
        name: managed group; reapplying replaces it {default: cuemol}
        edge: auto, none, edges, silhouette, thin, normal, thick
        edge_width: positive line width in angstroms (minimum one GPU pixel)
        edge_color: a PyMOL color name {default: black}
        transparency: keep, or 0 (opaque) to 1 (invisible) {default: keep}
        cache_mb: limit for each mesh and native CGO cache in MiB {default: 2048}
        filename: required output PNG path for png/ray
        width, height: output pixels; 0 preserves current dimensions

    NOTES

        All loaded states are prepared before playback. Use refresh after
        coordinate, topology, color, secondary-structure, or state edits.
        Richardson uses a paper-and-colored-pencil GPU shader. Native ray
        and transparent bodies approximate its stroke coverage as vertex tones.
        Opaque bodies and edges use GLSL. Transparent bodies use native CGO
        with baked lighting. Retained ray-only CGO also supports standard ray
        and png, ray=1 with approximate materials. Use the dedicated ray
        operation to add camera-dependent outlines and material samples.
        Reset restores native representations.
        Background settings are left to PyMOL and are never changed here.
        Use name=all to reset or refresh all managed views.
        Maps, labels, and unsupported representations remain native.
    """
    if _self is None:
        from pymol import cmd as _self
    if not _self.is_gui_thread() and callable(
        getattr(_self, "_call_in_gui_thread", None)
    ):
        parameters = locals().copy()
        return _self._call_in_gui_thread(lambda: cuemol_style(**parameters))
    from .controller import manager_for
    from .presets import PROFILES

    manager = manager_for(_self)
    manager.maintenance()
    style = str(style).strip().lower()
    try:
        if style == "list":
            print(" cuemol_style profiles: " + ", ".join(PROFILES))
            for entry in manager.entries.values():
                states = max(len(ds) for ds in entry.drawings.values())
                print(
                    f" {entry.name}: {entry.options['style']}, {states} states, prepared {entry.seconds:.2f} s, meshes {entry.nbytes / 1024**2:.1f} MiB, native CGO {entry.cgo_nbytes / 1024**2:.1f} MiB"
                )
            return tuple(PROFILES)
        if style == "reset":
            manager.reset(str(name))
            return
        if style == "refresh":
            manager.refresh(str(name))
            return
        if style in ("png", "ray"):
            from .export import image

            return image(manager, str(filename), width, height, style == "ray")
        entry = manager.apply(
            style,
            str(selection),
            None if representation is None else str(representation),
            str(color),
            str(quality),
            str(name),
            str(edge),
            None if edge_width in ("", None) else float(edge_width),
            edge_color,
            None if transparency in ("keep", "", None) else float(transparency),
            cache_mb=float(cache_mb),
        )
        if not int(quiet):
            states = max(len(ds) for ds in entry.drawings.values())
            print(
                f" cuemol_style: {name} ({style}), {states} states prepared in {entry.seconds:.2f} s; {entry.nbytes / 1024**2:.1f} MiB meshes, {entry.cgo_nbytes / 1024**2:.1f} MiB native CGO."
            )
        return entry
    except (ValueError, RuntimeError) as exc:
        from pymol import CmdException

        raise CmdException(str(exc)) from exc


def __init_plugin__(app=None):
    from pymol import cmd
    from pymol.shortcut import Shortcut

    from .controller import manager_for
    from .presets import COLORS, PROFILES, REPRESENTATIONS

    cmd.extend("cuemol_style", cuemol_style)
    cmd.auto_arg[0]["cuemol_style"] = [
        Shortcut([*PROFILES, "list", "refresh", "reset", "png", "ray"]),
        "style or operation",
        "",
    ]
    cmd.auto_arg[1]["cuemol_style"] = [cmd.selection_sc, "selection", ""]
    cmd.auto_arg[2]["cuemol_style"] = [Shortcut(REPRESENTATIONS), "representation", ""]
    cmd.auto_arg[3]["cuemol_style"] = [Shortcut(COLORS), "color mode", ""]
    manager_for(cmd)
