# cuemol_style_in_pymol

Standalone `cuemol_style` command for PyMOL: 26 CueMol-inspired molecular
geometry and material presets, interactive GPU rendering, picking, MD state
playback, transparent geometry, and standard or dedicated ray export.
PyMOL's source and standard commands are unchanged. CueMol is not required.

## Install

Use the reproducible environment in this repository:

```sh
pixi install
pixi run pymol
```

To install into an existing Python environment that already provides PyMOL 3.1:

```sh
uv pip install "git+https://github.com/th2ch-g/cuemol_style_in_pymol.git"
```

The Python package requires Python 3.10+, NumPy, SciPy, and PyOpenGL.
Interactive rendering additionally requires PyMOL 3.1 with Qt and a compatibility
OpenGL 2.1 / GLSL 1.20 context. PyMOL is supplied by conda/pixi, not by pip.

Register the command from PyMOL's Python console or add this to `.pymolrc.py`:

```python
from cuemol_style_in_pymol import __init_plugin__

__init_plugin__()
```

Importing the package alone does not apply a style or change PyMOL settings.
mdtbx's `pymol_plugins` integration registers the command automatically.

## Use

```text
cuemol_style richardson
cuemol_style richardson, representation=cpk
cuemol_style toon1
cuemol_style matte, representation=surface, transparency=0.4
ray 2400, 1800
png figure_ray.png
cuemol_style png, filename=figure.png, width=2400, height=1800
cuemol_style ray, filename=outlined_ray.png, width=2400, height=1800
cuemol_style refresh
cuemol_style reset
cuemol_style list
```

Default colors follow CueMol GUI's initial painting: khaki helices, SteelBlue
sheets, FloralWhite coils, and yellow nucleic geometry. Atomic representations
use DefaultCPKColoring, with carbon inheriting the molecular painting.
Helix outside faces keep their base color and inside faces use a lighter color.
The existing background is preserved. `color=keep` uses existing atom colors.

Standard `ray` and `png, ray=1` use retained native CGO. The dedicated ray
operation adds camera-dependent outline geometry and material samples.
The `richardson` GPU profile uses warm paper, three layers of irregular
colored-pencil strokes, bright unmarked highlights, and dark contour lines.
Standard ray, dedicated ray, and transparent CGO approximate the marks by their
average tone; they cannot run the per-pixel pencil shader.

See the [full guide](docs/pymol_cuemol.rst) for styles, states, selection,
session restoration, memory limits, and rendering differences.

## Gallery

Actual PyMOL GPU PNGs at `quality=high`, with CueMol default colors and a white
background. Protein views use [crambin, PDB 1CRN](https://www.rcsb.org/structure/1CRN).
DNA uses PyMOL's `fnab` builder; atomic views use its tryptophan fragment.
Run `cuemol_style <style>` for a geometry preset. Sticks use
`cuemol_style default, representation=sticks`.

| Ribbon | Round ribbon | Fancy ribbon |
| --- | --- | --- |
| ![Ribbon](docs/gallery/ribbon.png) | ![Round ribbon](docs/gallery/round_ribbon.png) | ![Fancy ribbon](docs/gallery/fancy_ribbon.png) |
| `ribbon` | `round_ribbon` | `fancy_ribbon` |
| ![Cartoon](docs/gallery/cartoon.png) | ![Round cartoon](docs/gallery/round_cartoon.png) | ![Tube](docs/gallery/tube.png) |
| `cartoon` | `round_cartoon` | `tube` |
| ![Nucleic](docs/gallery/nucleic.png) | ![Ball and stick](docs/gallery/ballstick.png) | ![Sticks](docs/gallery/sticks.png) |
| `nucleic` | `ballstick` | `representation=sticks` |
| ![CPK](docs/gallery/cpk.png) | ![Surface](docs/gallery/surface.png) | |
| `cpk` | `surface` | |

| Richardson ribbon | Richardson CPK |
| --- | --- |
| ![Richardson ribbon](docs/gallery/richardson.png) | ![Richardson CPK](docs/gallery/richardson_cpk.png) |
| `cuemol_style richardson` | `cuemol_style richardson, representation=cpk` |

The representation audit in the [full guide](docs/pymol_cuemol.rst) records the
matched dimensions and remaining differences. This is an independent renderer;
the images are not claimed to be pixel-identical to CueMol.

Regenerate the gallery in the repository's PyMOL Qt environment:

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/render_gallery.py
```

The published PNGs are versioned for README display. Downloaded coordinates and
temporary validation output stay in the ignored `.cache` directory.

## Develop and validate

```sh
pixi run test
pixi run check
pixi run -e docs docs
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --gui --benchmark --output .cache/validation
```

The standalone harness verifies real PyMOL ray export, all presets, restoration,
state changes, Qt picking, native/custom transparency, and GPU state preservation.
The optional benchmark uses 500 residues and 100 synthetic states.
Temporary images, environments, build products, and caches are ignored.
