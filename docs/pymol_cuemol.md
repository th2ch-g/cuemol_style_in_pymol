# CueMol-inspired PyMOL styles

[Package README](../README.md) | [日本語](ja/pymol_cuemol.md)

`cuemol_style_in_pymol` is a standalone package that adds molecular
geometry and live GPU materials to PyMOL. It requires PyMOL 3.1 with Qt,
NumPy, SciPy, and PyOpenGL, and a compatibility OpenGL 2.1 / GLSL 1.20
context. `pixi install` installs the development environment, including
PyMOL. CueMol and mdtbx are not required. PyMOL's source code and standard
commands are unchanged.

## Quick start

Run `pixi install` in this repository, or install the package into an
existing PyMOL Python environment with `uv pip install
"git+https://github.com/th2ch-g/cuemol_style_in_pymol.git"`. Register the
command from PyMOL's Python console or add these lines to `.pymolrc.py`:

```python
from cuemol_style_in_pymol import __init_plugin__

__init_plugin__()
```

Importing the package alone does not change PyMOL settings. The optional
mdtbx `pymol_plugins` integration registers the command automatically.
After loading a structure, run these commands in the PyMOL command line:

```text
cuemol_style richardson
cuemol_style toon1
cuemol_style toon2, selection=chain A
cuemol_style matte, representation=surface, transparency=0.4
cuemol_style list
help cuemol_style
```

The default named view is `cuemol`. Applying another style with the same
name replaces that view. `richardson` uses thin ribbons, sheet arrows,
lighter helix undersides, colored-pencil hatching, and black silhouette/crease lines.
Rotation and zoom update the GPU rendering immediately.

The `richardson` GPU shader follows CueMol 3's Richardson tone recipe:
warm paper, three irregular pencil layers at 55, -35, and 80 degrees,
pigment-colored strokes, and unmarked highlights. Tone depends on the
surface normal and viewing direction. Helix outside faces retain the base
pigment; inside faces are lighter. The live shader filters a fine pencil
lattice using CueMol's default 3x ink-sampling scale. Independent procedural
noise and the omission of Umbreon's occlusion, depth fog, and screen-space
stroke-edge pass make this an interactive approximation. The material and
outline audit below records the matched parameters and remaining differences.

## Styles and controls

Style names

| Group | Names |
| --- | --- |
| Protein geometry | `richardson`, `ribbon`, `round_ribbon`, `fancy_ribbon`, `cartoon`, `round_cartoon`, `tube` |
| Other molecular geometry | `nucleic`, `ballstick`, `cpk`, `surface` |
| Lighting and shading | `default`, `shadow`, `nolighting`, `matte`, `toon1`, `toon2` |
| Decorative materials | `diff_metal`, `spec_metal`, `metallic_chrome`, `metallic_copper`, `stone35`, `wood31`, `wood14scl2` |
| Outlines | `outline`, `silhouette` |

Geometry presets choose a representation. Material presets use
`representation=auto`: existing sticks, spheres, surface, cartoon, and
ribbon layers are retained as custom geometry. Atoms shown only as lines or
nonbonded points become protein ribbons, nucleic backbones/base-pair rods, or
ball-and-stick geometry. Hidden atoms stay hidden in this automatic mode.
Maps, labels, and other unsupported layers remain native.

Override the representation with `ribbon`, `cartoon`, `tube`,
`nucleic`, `ballstick`, `sticks`, `cpk`, or `surface`. `cartoon`
uses helix cylinders; `ribbon` follows the backbone. Quality is `low`,
`medium` (default), or `high`. Higher quality increases preparation time
and memory use.

`edge` accepts `auto`, `none`, `edges`, `silhouette`, `thin`,
`normal`, or `thick`. `edges` includes sharp visible creases.
`edge_width` is in angstroms, with a one-pixel minimum in GPU images;
`edge_color` accepts a PyMOL color. `transparency=keep` preserves the
source representation's opacity; a number from 0 to 1 overrides it.

```text
cuemol_style spec_metal, representation=cpk, edge=silhouette
cuemol_style richardson, edge_width=0.12, edge_color=black
cuemol_style wood31, representation=surface, quality=high
```

## Default colors

`color=cuemol` follows CueMol GUI's initial painting for newly loaded
molecules. Protein ribbons/cartoons/tubes use DefaultHSCPaint, nucleic
geometry uses DefaultNucl, and atomic/surface representations use
DefaultCPKColoring. These styles reference the molecular painting; carbon
inherits that painting in atomic/surface representations as well.

CueMol palette

| Feature | Color |
| --- | --- |
| Helix / sheet / coil | khaki `#F0E68C` / SteelBlue `#4682B4` / FloralWhite `#FFFAF0` |
| Nucleic backbone and bases | yellow `#FFFF00` |
| Carbon / nitrogen / oxygen | molecular painting / blue / red |
| Hydrogen / sulfur / phosphorus | cyan / green / yellow |

Other modes are `keep` (existing atom colors), `chain`, `ss`
(WoodyHSC secondary-structure colors), `rainbow`, and `element` (a
conventional element palette). Changing a PyMOL atom color while a style is
active takes effect after `refresh` and only when the color mode uses that
atom color. Source atom colors and secondary-structure assignments are
never changed by this command.

```text
cuemol_style richardson, color=keep
cuemol_style nucleic, color=chain
```

## Selection and state management

Click custom geometry to select its original atoms in `sele`. PyMOL's
mouse selection mode controls atom/residue/chain expansion; Shift adds to
the selection. Pink markers identify selected source atom positions.
Dragging retains native rotation/movement. Native opaque geometry blocks
clicks on custom geometry behind it. Picking requires the Qt GUI.

All loaded states are prepared before playback. Global state changes,
movie frame-to-state mappings, object state overrides, and `all_states`
are supported. Coordinate, bond, color, secondary-structure, transparency,
or state-count edits require `refresh`. To change the native representation
layers, reset the view, change the layers, and apply the style again. Source
visibility and object state settings are synchronized by a short Qt timer.

```text
cuemol_style refresh
cuemol_style reset
cuemol_style reset, name=all
```

`reset` restores the original native representation masks and removes
generated objects and GUI hooks. Applying, refreshing, resetting, and
restoring a style never change the background; use PyMOL's `bg_color` to
choose it. Source coordinates and bonds are untouched. Use distinct `name`
values for disjoint selections;
overlapping managed selections are rejected. Deleting a managed group or
source object triggers cleanup. Saved PyMOL sessions store recipes and
rebuild the views on load when the plugin is available.

`cache_mb=2048` separately limits prepared meshes, native CGO command
payloads, and estimated input snapshots. It is not a cap on total process
memory: each category has its own limit, and native allocations add overhead. The
GPU cache retains up to 256 MiB of vertex buffers and evicts old states;
evicted buffers are uploaded again from prepared CPU meshes. Large
trajectories, especially surfaces and transparent atomistic views, can be
expensive to prepare. `cuemol_style list` reports preparation time and
mesh/native CGO storage for each active view. Retaining ray geometry for all
states increases preparation time, memory use, and saved session size.

## Image export

```text
ray 2400, 1800
png figure_ray.png
png figure_ray.png, width=2400, height=1800, ray=1
cuemol_style png, filename=figure.png, width=2400, height=1800
cuemol_style ray, filename=figure_ray.png, width=2400, height=1800
```

Standard PyMOL `ray` and `png, ray=1` work directly, including in
headless mode. Opaque meshes have a retained CGO using PyMOL's ray-only
triangle opcode. These triangles do not draw in OpenGL or interfere with
picking. Transparent meshes already have native CGO and are not duplicated.
Every prepared state has matching ray geometry; global state and movie
frame changes take effect immediately. Source object visibility and
per-object state overrides use the same Qt synchronization as the live
view. Run `cuemol_style refresh` after changing those settings in headless
mode, or before an immediate ray command that cannot wait for the Qt timer.

The two `cuemol_style` export operations save the complete visible scene
and require a filename. `png` captures the GPU appearance and requires
the Qt GUI. `ray` builds temporary CGO geometry for the current state and
camera, adding explicit silhouette/crease lines and camera-dependent
material samples. It also works in headless PyMOL. Export restores
visibility, playback, and temporary settings even if rendering fails.
Selection markers are omitted from exports.

Opaque bodies use independent GPU shaders. Transparent bodies use native
CGO with lighting baked into vertex colors so that PyMOL can composite
them with native translucent objects. Their lighting is fixed in molecular
coordinates during rotation. Standard ray uses the same molecular-space
material samples and the user's PyMOL lighting/outline settings. The
dedicated ray operation instead uses camera-space samples and cylindrical
outline geometry. Both ray paths differ from the GPU shading and line
appearance. For Richardson, ray and transparent CGO use the average
colored-pencil coverage as vertex tones; individual hatching strokes are
available in opaque GPU rendering and `cuemol_style png`. Wood, stone, and
metal are procedural approximations; they do not reproduce CueMol's
POV-Ray textures exactly. `shadow` is a flat shading material, not a
scene-shadow generator.

Object transformation matrices, stereo/VR picking, editing atoms
by dragging, and headless interactive GPU rendering are outside the
supported interface; apply coordinate transforms to source atoms and
refresh when needed.

## Representation audit

Geometry defaults were checked against [CueMol 3173d8a](https://github.com/CueMol/cuemol2/tree/3173d8af62e211dd37b943ee53b3d6a632e6b5d7),
including `default_style.xml`, `TubeSection`, `RibbonRenderer`,
`Ribbon2Renderer`, `NARenderer`, and the atomic renderers. Richardson
tone parameters were checked against [Umbreon bf75c8a](https://github.com/CueMol/umbreon/tree/bf75c8adc05ed70a1344afbd718bcaab651c1070).
These are source-level checks; no claim of pixel equality is made.

- `ribbon` and `round_ribbon`: helix half-width 1.2, sheet half-width
  1.4, half-thickness 0.2, coil radius 0.35 angstrom. Axes use natural cubic
  splines with chord-length knots and 50 percent sheet-pivot smoothing.
  Sheet-arrow expansion is 1.8, with gamma 2.2 or 1.2 respectively.
- `fancy_ribbon` and `richardson`: helix half-width 1.3, circular rail
  radius 0.2, sharpness 0.3, sheet half-width 1.2, and coil radius 0.25.
  Helix backs and sheet side walls reduce HSV saturation by 0.4. Sheet
  arrows expand by 1.6 with gamma 1.0.
- `cartoon` and `round_cartoon`: penalized natural-spline helix axes
  use rho 3.0; cylinder radius is the mean pivot-to-axis distance plus
  0.2. Sheet half-width/thickness are 1.4/0.2, with smoothing rho 3.0 or
  1.0 respectively; coil radius is 0.2. Sheet arrows use expansion 1.8
  and gamma 1.0. Junction and endpoint constraints are approximations.
- `tube`: radius 0.35 with a natural cubic axis. `nucleic`: P-atom
  pivots, an elliptical backbone with half-axes 1.25/0.5, and base-pair
  rods of radius 0.5. Base pairs are inferred from compatible in-plane
  hydrogen-bond contacts. Pair assignment and modified-base support can
  differ from CueMol's residue topology and base-pair metadata.
- `ballstick`: every atom radius 0.3 and bond radius 0.2. `sticks`:
  atom and bond radius 0.2. Bond colors split sharply at the midpoint.
  `cpk`: H/C/N/O/S/P radii 1.2/1.7/1.55/1.52/1.8/1.8, other elements
  1.7. Mesh tessellation depends on this plugin's quality setting.
- `surface`: solvent-excluded surface with a 1.4-angstrom probe and
  the same element radii. PyMOL's surface mesher differs from CueMol's
  EDTSurf/MeshMS implementation, so triangulation and fine details differ.

Spline frames, chain-break detection, section transitions, caps, and
surface ownership use this plugin's implementation. Matching default
dimensions does not make every molecular representation identical.
Maps and labels retain their existing native PyMOL representation.

## Material and outline audit

The OpenGL material coefficients match CueMol's `default_style.xml`:

| Material | Ambient | Diffuse | Specular | Shininess |
| --- | ---: | ---: | ---: | ---: |
| `default` | 0.2 | 0.8 | 0.0 | 32.0 |
| `shadow` | 0.75 | 0.0 | 0.0 | 0.0 |
| `nolighting` | 1.0 | 0.0 | 0.0 | 0.0 |
| `matte` | 0.3 | 0.6 | 0.0 | 32.0 |
| `toon1` | 0.0 | 0.85 | 0.0 | 0.0 |
| `toon2` | 0.0 | 0.85 | 0.0 | 32.0 |
| `diff_metal`, `spec_metal` | 0.2 | 0.5 | 0.7 | 76.8 |

GPU and native material samples share these coefficients. The camera-space
light direction is `(1, 1, 1.5)`, following CueMol's OpenGL lighting source.
`toon1` and `toon2` consequently have the same diffuse appearance, as do
`diff_metal` and `spec_metal`. CueMol's POV-Ray definitions distinguish those
pairs using different finishes; this renderer does not evaluate POV-Ray
`brilliance`, `phong`, or `F_MetalA`/`F_MetalD`. CueMol's backend lighting,
occlusion, shadows, and display transfer can still produce different images.

`metallic_chrome`, `metallic_copper`, `stone35`, `wood31`, and `wood14scl2`
approximate the corresponding POV-Ray texture names. Their procedural texture
patterns and reflection bands are this plugin's own implementation; those
settings and pixels are not identical to CueMol's POV-Ray texture library.

The `thin`, `normal`, and `thick` outline widths are 0.03, 0.06, and
0.15 angstrom, matching CueMol's named EgLine styles. `outline` enables
silhouettes and sharp creases; `silhouette` enables silhouettes alone.
The combined `richardson`, `toon1`, and `toon2` profiles select normal black
edges. GPU edge extraction and its one-pixel minimum, and cylindrical ray
outlines, differ from CueMol's image-space stroke rendering. CueMol's raw
renderer default of 0.01 angstrom is separate from its named normal edge style.

Richardson uses the same khaki and SteelBlue pigments, paper `(0.941, 0.925,
0.867)`, layer angles `55/-35/80`, thresholds `0.92/0.62/0.34`, width fade `10`,
ink scales `1/0.74/0.38`, dark-pressure floor `0.4`, and minimum contrast `0.15`.
The configured spacing/width are `0.5/0.45` output pixels. With CueMol's default
3x supersampling and its two-device-pixel minimum, the effective pitch is
`2/3` output pixel. The live shader now uses that pitch and width, a nine-sample
stroke filter, length/gap `50/5`, width jitter `0.45`, length jitter `0.5`,
taper `0.35`, angle jitter `5` degrees, and paper tooth `0.15` at scale `3`.
The shader's noise is independent; it shares the slow stroke envelope across
the nine samples and averages each layer before multiplying the layers.
Umbreon evaluates and multiplies them at every supersample. Occlusion and
depth fog are also omitted. Highlights therefore remain bare paper, and the
result is not a pixel-identical Umbreon render. Native ray/transparency use
a fitted average coverage, about 0.256 per fully active pencil layer.

## Validation

The standalone harness exercises every style in real PyMOL, restoration
after failures, multiple states, session reload, Qt picking, rotation, and
native/custom transparency. GUI mode also writes GPU and ray images for
visual review. Run Python through the repository's pixi interpreter:

```console
$ uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --output .cache/cuemol-headless
$ uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --gui --benchmark \
    --output .cache/cuemol-gui
```

The benchmark uses 500 residues and 100 synthetic states at a 1280 by 720
viewport with medium-quality ribbons. It reports preparation time, CPU/GPU
mesh storage, rotation speed, explicit state-switch speed, actual movie
draw rate, and one-state surface preparation. The targets after preparation
are 30 FPS rotation and 15 FPS playback; results depend on the GPU and input
geometry. Generated images and reports under `.cache` are ignored by Git.

The geometry and names are inspired by the
[CueMol ribbon renderer](https://cuemol.github.io/cuemol2_docs/cuemol2/RibbonRenderer/)
and CueMol style definitions. This plugin implements its own mesh generation
and rendering and has no runtime dependency on the CueMol source tree.
