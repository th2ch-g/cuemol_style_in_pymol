"""Render README gallery images with a separate real PyMOL Qt process."""

import argparse
import sys
from pathlib import Path
from time import monotonic, sleep
from urllib.request import urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/gallery"))
    parser.add_argument("--cache", type=Path, default=Path(".cache/gallery"))
    parser.add_argument(
        "--only", nargs="+", help="Render only the named gallery entries"
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    structure = args.cache / "1crn.pdb"
    if not structure.exists():
        with urlopen(
            "https://files.rcsb.org/download/1CRN.pdb", timeout=30
        ) as response:
            structure.write_bytes(response.read())
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    import pymol

    from cuemol_style_in_pymol import cuemol_style
    from cuemol_style_in_pymol.presets import PROFILES

    pymol.invocation.options.show_splash = 0
    from pmg_qt.pymol_qt_gui import PyMOLApplication, PyMOLQtGUI

    app = PyMOLApplication(["cuemol-gallery"])
    window = PyMOLQtGUI()
    widget, cmd = window.pymolwidget, window.pymolwidget.cmd

    def in_context(func):
        with widget:
            return func()

    def pump(seconds=0.1):
        end = monotonic() + seconds
        while monotonic() < end:
            app.processEvents()
            sleep(0.002)

    cmd._call_with_opengl_context = in_context
    window.resize(1100, 900)
    window.show()
    pump(0.3)
    try:
        cmd.set("internal_gui", 0)
        cmd.set("internal_feedback", 0)
        for setting in ("seq_view", "movie_panel", "internal_prompt"):
            cmd.set(setting, 0)
        ratio = widget.devicePixelRatioF()
        widget.setFixedSize(round(1200 / ratio), round(900 / ratio))
        window.adjustSize()
        pump(0.3)
        with widget:
            widget.resizeGL(widget.width(), widget.height())
        cmd.set("orthoscopic", 1)
        cmd.set("antialias", 2)
        cmd.bg_color("white")
        groups = [
            (
                "protein",
                [
                    "ribbon",
                    "round_ribbon",
                    "fancy_ribbon",
                    "cartoon",
                    "round_cartoon",
                    "tube",
                    "surface",
                    "richardson",
                ],
            ),
            ("dna", ["nucleic"]),
            ("atoms", ["ballstick", "sticks", "cpk", "richardson_cpk"]),
        ]
        covered = {name for _, names in groups for name in names}
        groups[0][1].extend(name for name in PROFILES if name not in covered)
        if args.only:
            groups = [
                (sample, [name for name in names if name in args.only])
                for sample, names in groups
            ]
            unknown = set(args.only) - {name for _, names in groups for name in names}
            if unknown:
                parser.error(f"Unknown gallery entries: {sorted(unknown)}")
        for sample, profiles in groups:
            if not profiles:
                continue
            cmd.delete("all")
            if sample == "protein":
                cmd.load(str(structure), sample)
                cmd.remove("not polymer.protein")
                cmd.dss(sample)
            elif sample == "dna":
                cmd.fnab("ATGCGCAT", name=sample)
            else:
                cmd.fragment("trp", sample)
            cmd.show_as("lines", sample)
            cmd.orient(sample)
            if sample == "protein":
                cmd.turn("y", 20)
                cmd.turn("z", -15)
            pump()
            cmd.zoom(sample, 3, complete=1)
            cmd.clip("slab", 150)
            view = cmd.get_view()
            for profile in profiles:
                cmd.set_view(view)
                print("Rendering", profile, flush=True)
                style = (
                    "richardson"
                    if profile == "richardson_cpk"
                    else "default"
                    if profile == "sticks"
                    else profile
                )
                representation = (
                    "cpk"
                    if profile == "richardson_cpk"
                    else "sticks"
                    if profile == "sticks"
                    else None
                )
                cuemol_style(
                    style,
                    representation=representation,
                    quality="high",
                    quiet=1,
                    _self=cmd,
                )
                pump()
                cuemol_style(
                    "png",
                    filename=str(args.output / f"{profile}.png"),
                    width=1200,
                    height=900,
                    _self=cmd,
                )
                cuemol_style("reset", _self=cmd)
    finally:
        cuemol_style("reset", _self=cmd)
        cmd.delete("all")
        window.hide()


if __name__ == "__main__":
    main()
