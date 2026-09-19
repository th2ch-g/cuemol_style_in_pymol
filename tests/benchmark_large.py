"""Measure preparation using a local structure in a fresh PyMOL process."""

import argparse
import json
import resource
import sys
from pathlib import Path
from time import perf_counter

import pymol2

from cuemol_style_in_pymol import cuemol_style


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structure", type=Path)
    parser.add_argument("--selection", default="all")
    parser.add_argument("--style")
    parser.add_argument(
        "--output", type=Path, default=Path(".cache/large-structure.json")
    )
    args = parser.parse_args()
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.load(str(args.structure), "structure", quiet=1)
        result = {
            "atoms": cmd.count_atoms(args.selection),
            "states": cmd.count_states(args.selection),
            "pymol": cmd.get_version()[0],
        }
        options = {"selection": args.selection, "quiet": 1, "_self": cmd}
        if args.style:
            options["style"] = args.style
        started = perf_counter()
        entry = cuemol_style(**options)
        result["prepare_seconds"] = perf_counter() - started
        if hasattr(entry, "nbytes"):
            result["geometry_bytes"] = entry.nbytes
            result["native_bytes"] = getattr(
                entry, "cgo_nbytes", getattr(entry, "native_bytes", 0)
            )
            result["compact_spheres"] = sum(
                len(part.spheres)
                for ds in entry.drawings.values()
                for d in ds
                for part in d.native
            )
            result["compact_cylinders"] = sum(
                len(part.cylinders)
                for ds in entry.drawings.values()
                for d in ds
                for part in d.native
            )
        result["peak_process_bytes"] = resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
        cuemol_style("reset", name="all", quiet=1, _self=cmd)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
