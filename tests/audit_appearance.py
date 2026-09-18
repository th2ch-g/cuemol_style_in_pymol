"""Render portable, camera-aligned PyMOL/CueMol comparison cases."""

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
from check_cuemol_style import qt_session
from scipy.spatial.transform import Rotation

from cuemol_style_in_pymol import cuemol_style, geometry
from cuemol_style_in_pymol.presets import PROFILES, resolve
from cuemol_style_in_pymol.source import read

REFERENCE_REVISION = "3173d8af62e211dd37b943ee53b3d6a632e6b5d7"


def secondary_records(atoms):
    residues = [a for a in atoms if a.name == "CA"]
    records = []
    start = 0
    while start < len(residues):
        end = start + 1
        first = residues[start]
        while (
            end < len(residues)
            and residues[end].ss == first.ss
            and residues[end].chain == first.chain
        ):
            end += 1
        if first.ss in ("H", "S"):
            last = residues[end - 1]
            line = list(" " * 80)
            if first.ss == "H":
                line[:6] = "HELIX "
                line[19], line[31] = first.chain or " ", last.chain or " "
                line[21:25], line[33:37] = (
                    f"{int(first.resi):4d}",
                    f"{int(last.resi):4d}",
                )
                line[38:40] = " 1"
            else:
                line[:6] = "SHEET "
                line[21], line[32] = first.chain or " ", last.chain or " "
                line[22:26], line[33:37] = (
                    f"{int(first.resi):4d}",
                    f"{int(last.resi):4d}",
                )
            records.append("".join(line))
        start = end
    return "\n".join(records) + "\n"


def reference_style(profile, representation):
    if representation == "ribbon":
        return "ribbon", {"fancy": "Fancy1Ribbon", "ellipse": "RoundRibbon"}.get(
            profile.section, "DefaultRibbon"
        )
    if representation == "cartoon":
        return (
            "cartoon",
            "RoundCartoon" if profile.section == "ellipse" else "DefaultCartoon",
        )
    return {
        "tube": ("tube", ""),
        "nucleic": ("nucl", "DefaultNucl"),
        "ballstick": ("ballstick", "DefaultBallStick"),
        "sticks": ("ballstick", "StickBallStick"),
        "cpk": ("cpk", "DefaultCPK"),
        "surface": ("dsurface", ""),
    }[representation]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True, type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path(".cache/appearance-audit/aligned")
    )
    parser.add_argument("--only", nargs="+", default=[*PROFILES, "sticks"])
    parser.add_argument("--reference-module", type=Path)
    parser.add_argument("--reference-config", type=Path)
    parser.add_argument("--angle", type=float, default=0)
    parser.add_argument("--zoom", type=float, default=1)
    parser.add_argument("--perspective", action="store_true")
    output_mode = parser.add_mutually_exclusive_group()
    output_mode.add_argument("--ray", action="store_true")
    output_mode.add_argument("--live", action="store_true")
    parser.add_argument("--transparency", type=float, default=0)
    parser.add_argument("--width", type=int, default=640)
    parser.add_argument("--height", type=int, default=480)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with qt_session() as (cmd, widget, pump):
        cmd.set("internal_gui", 0)
        cmd.set("opaque_background", 1)
        cmd.set("ray_opaque_background", 1)
        for setting in (
            "internal_feedback",
            "seq_view",
            "movie_panel",
            "internal_prompt",
        ):
            cmd.set(setting, 0)
        ratio = widget.devicePixelRatioF()
        widget.setFixedSize(round(args.width / ratio), round(args.height / ratio))
        pump(0.2)
        with widget:
            widget.resizeGL(widget.width(), widget.height())
        for name in args.only:
            cuemol_style("reset", quiet=1, _self=cmd)
            cmd.delete("all")
            profile = (
                resolve("default", "sticks") if name == "sticks" else resolve(name)
            )
            representation = profile.representation
            if name == "nucleic":
                cmd.fnab("ATGCGCAT", name="molecule")
            elif representation in ("cpk", "ballstick", "sticks"):
                cmd.fragment("trp", "molecule")
            else:
                cmd.load(str(args.structure), "molecule")
                cmd.dss("molecule")
            cmd.remove("hydro or solvent")
            cmd.show_as("lines")
            cmd.orient("molecule")
            cmd.turn("y", 20 + args.angle)
            cmd.turn("z", -15)
            cmd.zoom("molecule", 3)
            if args.zoom != 1:
                view = list(cmd.get_view())
                view[11] /= args.zoom
                cmd.set_view(view)
            cmd.clip("slab", 150)
            cmd.set("orthoscopic", not args.perspective)
            background = (
                [240 / 255, 236 / 255, 221 / 255] if name == "richardson" else [1, 1, 1]
            )
            cmd.set_color("audit_background", background)
            cmd.bg_color("audit_background")
            before, _ = read(cmd, "molecule", 2**30)
            frame = before["molecule"][0]
            entry = cuemol_style(
                "default" if name == "sticks" else name,
                representation=representation,
                transparency=args.transparency,
                quality="medium",
                quiet=1,
                _self=cmd,
            )
            pump()
            image_path = str(args.output / f"{name}-pymol.png")
            if args.live:
                with widget:
                    widget.paintGL()
                if not widget.grabFramebuffer().save(image_path):
                    raise RuntimeError("Could not save the live framebuffer")
            else:
                cuemol_style(
                    "ray" if args.ray else "png",
                    filename=image_path,
                    width=args.width,
                    height=args.height,
                    quiet=1,
                    _self=cmd,
                )
            drawing = entry.drawings["molecule"][0]
            if drawing.error or drawing.matrices is None:
                raise RuntimeError(drawing.error or "The drawing callback did not run")
            modelview, projection, _viewport = drawing.matrices
            if args.ray:
                from cuemol_style_in_pymol.sampling import camera

                modelview, projection, _viewport = camera(
                    cmd, args.width, args.height, ray=True
                )
            coordinates = frame.coords
            view = np.asarray(cmd.get_view())
            distance = float(-view[11])
            colors = geometry.atom_colors(frame.atoms, "cuemol", representation)
            pdb_name = f"{name}.pdb"
            cmd.save(str(args.output / pdb_name), "molecule")
            pdb_path = args.output / pdb_name
            pdb_path.write_text(secondary_records(frame.atoms) + pdb_path.read_text())
            renderer, styles = reference_style(profile, representation)
            manifest = {
                "reference_revision": REFERENCE_REVISION,
                "structure": pdb_name,
                "output": f"{name}-cuemol.png",
                "renderer": renderer,
                "styles": styles,
                "properties": {
                    "material": profile.material
                    if profile.material != "richardson"
                    else "default",
                    "egtype": profile.edges,
                    "eglinew": profile.edge_width,
                    "alpha": 1 - args.transparency,
                },
                "background": background,
                "rotation": Rotation.from_matrix(modelview[:3, :3].T)
                .as_quat()
                .tolist(),
                "center": view[12:15].tolist(),
                "zoom": float(2 / projection[1, 1])
                * (distance if args.perspective else 1),
                "distance": distance,
                "slab": 150,
                "width": args.width,
                "height": args.height,
                "perspective": args.perspective,
                "hatching": name == "richardson",
                "transparent_background": False,
                "pymol_export": "GPU live"
                if args.live
                else "ray"
                if args.ray
                else "GPU PNG",
                "atoms": [
                    {
                        "name": a.name,
                        "resi": a.resi,
                        "chain": a.chain,
                        "ss": a.ss,
                        "position": p.tolist(),
                        "color": c.tolist(),
                    }
                    for a, p, c in zip(frame.atoms, coordinates, colors)
                ],
            }
            manifest_path = args.output / f"{name}.json"
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
            if args.reference_module and args.reference_config:
                result = subprocess.run(
                    [
                        "node",
                        str(Path(__file__).with_name("appearance_reference.cjs")),
                        str(args.reference_module),
                        str(args.reference_config),
                        str(manifest_path),
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                (args.output / f"{name}-reference.log").write_text(
                    result.stdout + result.stderr
                )
                result.check_returncode()
            print(f"Rendered {name}", flush=True)


if __name__ == "__main__":
    main()
