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

The Python package requires Python 3.10+, NumPy, SciPy, PyOpenGL, and Pillow.
Building from source also requires a C++17 compiler; wheels contain the
standalone EDTSurf, contour, and pencil-sampling extension.
Interactive rendering additionally requires PyMOL 3.1 with Qt and a compatibility
OpenGL 2.1 / GLSL 1.20 context. PyMOL is supplied by conda/pixi, not by pip.

Register the command from PyMOL's Python console or add this to `.pymolrc.py`:

```python
from cuemol_style_in_pymol import __init_plugin__

__init_plugin__()
```

Importing the package alone does not apply a style or change PyMOL settings.
mdtbx's `pymol_plugins` integration registers the command automatically.

## Gallery

All 26 styles, plus explicit sticks and Richardson CPK views. Every PNG is
1200 x 900 at `quality=high`, with CueMol default colors and a white background.
Protein views share the same camera and use
[crambin, PDB 1CRN](https://www.rcsb.org/structure/1CRN). DNA uses PyMOL's
`fnab` builder; atomic views use its tryptophan fragment.

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/ribbon.png" alt="ribbon" width="240"><br><code>ribbon</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/round_ribbon.png" alt="round_ribbon" width="240"><br><code>round_ribbon</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/fancy_ribbon.png" alt="fancy_ribbon" width="240"><br><code>fancy_ribbon</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/cartoon.png" alt="cartoon" width="240"><br><code>cartoon</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/round_cartoon.png" alt="round_cartoon" width="240"><br><code>round_cartoon</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/tube.png" alt="tube" width="240"><br><code>tube</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/nucleic.png" alt="nucleic" width="240"><br><code>nucleic</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/ballstick.png" alt="ballstick" width="240"><br><code>ballstick</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/sticks.png" alt="sticks" width="240"><br><code>default</code><br><code>representation=sticks</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/cpk.png" alt="cpk" width="240"><br><code>cpk</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/surface.png" alt="surface" width="240"><br><code>surface</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/richardson.png" alt="richardson" width="240"><br><code>richardson</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/richardson_cpk.png" alt="richardson_cpk" width="240"><br><code>richardson</code><br><code>representation=cpk</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/default.png" alt="default" width="240"><br><code>default</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/shadow.png" alt="shadow" width="240"><br><code>shadow</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/nolighting.png" alt="nolighting" width="240"><br><code>nolighting</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/matte.png" alt="matte" width="240"><br><code>matte</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/toon1.png" alt="toon1" width="240"><br><code>toon1</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/toon2.png" alt="toon2" width="240"><br><code>toon2</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/diff_metal.png" alt="diff_metal" width="240"><br><code>diff_metal</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/spec_metal.png" alt="spec_metal" width="240"><br><code>spec_metal</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/metallic_chrome.png" alt="metallic_chrome" width="240"><br><code>metallic_chrome</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/metallic_copper.png" alt="metallic_copper" width="240"><br><code>metallic_copper</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/stone35.png" alt="stone35" width="240"><br><code>stone35</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/wood31.png" alt="wood31" width="240"><br><code>wood31</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/wood14scl2.png" alt="wood14scl2" width="240"><br><code>wood14scl2</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/outline.png" alt="outline" width="240"><br><code>outline</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/silhouette.png" alt="silhouette" width="240"><br><code>silhouette</code></td>
    <td width="33%"></td>
    <td width="33%"></td>
  </tr>
</table>

Run `cuemol_style <style>` for a named preset. The
[full guide](docs/pymol_cuemol.md#representation-audit) records matched geometry
dimensions, material coefficients, and remaining rendering differences.
The comparison target is CueMol's current Umbreon direct renderer. The toon
and metal pairs have distinct finishes. Wood and stone names use its PBR
settings; legacy POV-Ray procedural textures belong to a different backend.

Regenerate the gallery in the repository's PyMOL Qt environment:

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/render_gallery.py
```

The published PNGs are versioned for README display. Downloaded coordinates
and temporary validation output stay in the ignored `.cache` directory.

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

Material and outline styles default to `ribbon`, keeping helices spiral-shaped
even when the source uses PyMOL's cartoon display. Geometry presets such as
`cpk`, `surface`, and `cartoon` select their named representation. Use
`representation=auto` explicitly to inherit the source layers with a material
or outline style.

Default colors follow CueMol GUI's initial painting: khaki helices, SteelBlue
sheets, FloralWhite coils, and yellow nucleic geometry. Atomic representations
use DefaultCPKColoring, with carbon inheriting the molecular painting.
Fancy/Richardson helices retain their pigment on the outside and rounded rails;
the flat underside is lighter.
The existing background is preserved. `color=keep` uses existing atom colors.

Standard `ray` and `png, ray=1` use retained native CGO. The dedicated ray
operation adds camera-dependent outline geometry and material samples.
The `richardson` GPU profile uses warm paper, three layers of irregular
colored-pencil strokes, bright unmarked highlights, and dark contour lines.
Opaque display uses a tiled 3x render pass. Screen-space depth and normal
continuity define joined contours; internal sheet triangulation is never drawn
as an outline. Transparent CGO and dedicated ray share visible-sample colors and
contours on the same 3x grid. Dedicated exports combine completed transparent
renderer passes in display space, with alpha coverage preserved. These paths favor
appearance over frame rate and can require substantial preparation time and
memory. Standard ray retains an averaged pencil fallback. Dedicated ray
temporarily neutralizes PyMOL lighting
for the whole exported scene; use GPU PNG when native objects must retain their
existing lighting.

See the [full guide](docs/pymol_cuemol.md) for styles, states, selection,
session restoration, memory limits, and rendering differences.
The [Japanese guide](docs/ja/pymol_cuemol.md) covers the same interface.

## Develop and validate

```sh
pixi run test
pixi run check
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --gui --benchmark --output .cache/validation
```

The standalone harness verifies real PyMOL ray export, all presets, restoration,
state changes, Qt picking, native/custom transparency, and GPU state preservation.
The optional benchmark uses 500 residues and 100 synthetic states.
The [appearance audit](docs/pymol_cuemol.md#validation) compares real CueMol and
PyMOL exports with identical atoms, colors, secondary structure, and cameras.
Temporary images, environments, build products, and caches are ignored.
