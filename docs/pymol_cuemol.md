# CueMol-inspired PyMOL styles

[Package README](../README.md) | [日本語](ja/pymol_cuemol.md)

`cuemol_style_in_pymol` is a standalone package that adds molecular
geometry and live GPU materials to PyMOL. It requires PyMOL 3.1 with Qt,
NumPy, SciPy, PyOpenGL, and Pillow, and a compatibility OpenGL 2.1 / GLSL 1.20
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
cuemol_style
cuemol_style richardson
cuemol_style toon1
cuemol_style toon2, selection=chain A
cuemol_style matte, representation=surface, transparency=0.4
cuemol_style list
help cuemol_style
```

With no arguments, `cuemol_style` applies the `ribbon` preset.
The default named view is `cuemol`. Applying another style with the same
name replaces that view. `richardson` uses thin ribbons, sheet arrows,
lighter helix undersides, colored-pencil hatching, and black contour lines.
Rotation and zoom update the GPU rendering immediately.

The `richardson` pass uses warm paper, three pencil layers at 55, -35, and
80 degrees, pigment-colored strokes, and unmarked highlights. Flat helix
undersides are lighter; rounded rails retain their pigment. Opaque interactive
display evaluates pencil strokes and depth/normal contours entirely on the GPU.
Tiled 3x sampling retains fine strokes and smooth silhouettes, avoiding CPU
rasterization and framebuffer downloads while rotating or playing states.
Plain materials without outlines draw directly.

Dedicated PNG and ray keep the precision path: normals and
pigments are rasterized at 3x, then a native sampler evaluates the reference
hash/noise, stroke envelopes, and per-sample layer multiplication. Tone includes
depth fog. Bounded render tiles prevent oversized framebuffer allocations.
Contours are extracted from the visible depth/normal image, traced across pixel
boundaries, joined at occlusion junctions, smoothed, and rasterized as round
bands. Internal triangulation edges do not participate in contour extraction.

Molecular selections use original source objects. Native copies owned by
`chimerax_style` are excluded, including after a session reload, so applying
or refreshing a style never treats those display copies as input structures.
Ordinary source objects whose names start with an underscore remain supported.
Each command manages its own views; reset the previous command before switching
renderers on the same atoms to replace its display and restore its settings.

## Styles and controls

Style names

| Group | Names |
| --- | --- |
| Protein geometry | `richardson`, `ribbon`, `round_ribbon`, `fancy_ribbon`, `cartoon`, `round_cartoon`, `tube` |
| Other molecular geometry | `nucleic`, `ballstick`, `cpk`, `surface` |
| Lighting and shading | `default`, `shadow`, `nolighting`, `matte`, `toon1`, `toon2` |
| Decorative materials | `diff_metal`, `spec_metal`, `metallic_chrome`, `metallic_copper`, `stone35`, `wood31`, `wood14scl2` |
| Outlines | `outline`, `silhouette` |

Geometry presets choose their named representation. Material and outline
presets default to `ribbon`, so helices remain spiral-shaped even when the
source is shown with PyMOL's cartoon representation. For example,
`cuemol_style toon1` uses ribbons; `cuemol_style cartoon` explicitly uses helix
cylinders. An explicit representation always overrides the default.

With a material or outline style, `representation=auto` retains existing
sticks, spheres, surface, cartoon, and ribbon layers as custom geometry.
Atoms shown only as lines or nonbonded points become protein ribbons, nucleic
backbones/base-pair rods, or ball-and-stick geometry. Hidden atoms stay hidden
in this explicit automatic mode. Maps, labels, and other unsupported layers
remain native.

Override the representation with `ribbon`, `cartoon`, `tube`,
`nucleic`, `ballstick`, `sticks`, `cpk`, or `surface`. `cartoon`
uses helix cylinders; `ribbon` follows the backbone. Quality is `low`,
`medium` (default), or `high`. Higher quality increases preparation time
and memory use.

`edge` accepts `auto`, `none`, `edges`, `silhouette`, `thin`,
`normal`, or `thick`. `edges` includes self-occlusion contours; `silhouette`
keeps outer contours. Crease lines are disabled, matching the current
reference exporter's default crease limit.
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

All loaded states' base meshes are prepared before playback. Transparent
material samples are rebuilt outside drawing callbacks when the camera,
viewport, fog settings, or active state changes. Global state changes,
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

The interactive framebuffer uses 24 bytes per supersample in a reusable 3x
tile: about 54.4 MiB for a 510-pixel tile with a two-pixel overlap. Wider outlines
increase the overlap; allocations are capped at 256 MiB and the driver's texture
limit. The separate export render tile uses at most about 126 MiB of GPU attachments plus
temporary CPU arrays. Contours use a full-view 3x depth/normal image to keep
chain connectivity independent of tile boundaries. Native camera samples can
be larger than the base mesh: each covered 1/3-pixel sample uses two triangles. Inactive states
return to coarse fallback CGO instead of accumulating dense samples. Exceeding
`cache_mb` reports an error without silently reducing quality.

## Image export

```text
ray 2400, 1800
png figure_ray.png
cuemol_style png, filename=figure.png, width=2400, height=1800
cuemol_style ray, filename=figure_ray.png, width=640, height=480
```

Standard PyMOL `ray` and `png, ray=1` use retained native geometry for every
prepared state, including in headless mode. Opaque meshes use the ray-only
triangle opcode, which does not draw in OpenGL or interfere with picking.
Transparent geometry already has native CGO and is not duplicated. Standard
ray uses coarse molecular-space material samples, an averaged Richardson tone,
and the user's PyMOL lighting. Refresh after source visibility or object-state
changes in headless mode, or before a ray command that cannot wait for the Qt
maintenance timer.

The dedicated `png` operation requires Qt, omits selection markers, preserves
the current background, and renders opaque bodies in bounded tiles at 3x.
Transparent renderer groups use complete, opaque GPU layer passes with prepared
color/depth textures; bodies and contour ink are sampled together. This avoids
native CGO's extra lighting and fog on exported ink. The passes use the same 3x
grid and are downsampled before display-space group blending. Transfer buffers,
pixel storage, framebuffers, viewport, matrices, shader programs, and GL
attributes are restored.

The dedicated `ray` operation builds camera-dependent pixel CGO, including
individual pencil strokes and the same joined screen contours as GPU output.
Each covered sample has one front-facing quad at its visible depth. Hidden
triangles and separate contour cylinders cannot accumulate extra opacity.
Ray uses the same 3x grid and group compositing, including in headless mode.
To avoid shading the samples twice, it temporarily sets
neutral lighting and disables ray shadows/fog for the **whole exported scene**.
Unmanaged objects in that image therefore also use neutral lighting. All
settings, visibility, and playback are restored afterward, including on failure.
Use GPU PNG when native objects must retain their existing lighting.

Live transparent bodies use native CGO for integration with standard PyMOL
objects. Pieces sharing an opacity are sampled as one visible surface, so
internal overlaps and contour ink contribute opacity once. Colors and fog
compensation follow the camera, while source transparency settings are preserved.
The live native pass retains PyMOL's inter-object sorting and blend rules.

Dedicated PNG/ray exports use Umbreon's renderer-group formula instead:
`(1 - sum(alpha_i)) * B + sum(alpha_i * L_i)`, where `B` excludes all transparent
groups and `L_i` includes the background scene plus group `i` rendered opaque.
RGB is accumulated in sRGB-encoded space after shading, contours, and 3x
downsampling; coverage alpha is accumulated linearly. Background weights can
be negative when the group opacity sum exceeds one. Unmanaged geometry is
present in each pass, so occlusion against objects behind transparent groups
is included. Temporary objects and visibility are restored on failure.

Object transformation matrices, stereo/VR picking, editing atoms by dragging,
and headless interactive GPU rendering are outside the supported interface.
Apply coordinate transforms to source atoms and refresh when needed.

## Representation audit

The reference definitions are [CueMol 3173d8a](https://github.com/CueMol/cuemol2/tree/3173d8af62e211dd37b943ee53b3d6a632e6b5d7)
and [Umbreon bf75c8a](https://github.com/CueMol/umbreon/tree/bf75c8adc05ed70a1344afbd718bcaab651c1070).
Both numerical definitions and actual exports are compared. The target is
the current Umbreon direct renderer, with GI, AO, and cast shadows disabled.
Legacy OpenGL, POV-Ray procedural materials, and GI rendering are different
reference modes; pixel equality across these backends is not claimed.

| Geometry | Reference dimensions and construction |
| --- | --- |
| `ribbon`, `round_ribbon` | Helix half-width 1.2, sheet half-width 1.4, half-thickness 0.2, coil radius 0.35 A. Natural chord-length splines, 50% sheet-pivot smoothing, transported frames, slope normals, gamma 2.2 junctions. Arrow expansion 1.8; arrow gamma 2.2/1.2. |
| `fancy_ribbon`, `richardson` | Helix half-width 1.3, rail radius 0.2, sharpness 0.3, sheet half-width 1.2, coil radius 0.25 A. HSV saturation reduction 0.4 only on flat helix backs and sheet sides. Arrow expansion 1.6, gamma 1.0. |
| `cartoon`, `round_cartoon` | Penalized natural splines including flanking residues. Helix rho 3.0, radius = mean pivot-axis distance + 0.2 A. Sheet half-width/thickness 1.4/0.2 A, rho 3.0/1.0, facing-vector rho 5.0. Coil radius 0.2 A, rho -1/-2, weighted anchors and sheet-end derivative support. |
| `tube`, `nucleic` | Tube radius 0.35 A. Nucleic P-atom backbone half-axes 1.25/0.5 A, base-pair rods radius 0.5 A. Spline terminal caps use five hemispherical rings. Compatible in-plane hydrogen bonds determine base pairs. |
| `ballstick`, `sticks`, `cpk` | Ball/stick radii 0.3/0.2 A; sticks use 0.2/0.2 A. CPK H/C/N/O/S/P radii 1.2/1.7/1.55/1.52/1.8/1.8 A, other elements 1.7 A. Bond colors split at the midpoint. Dense spheres/cylinders approximate analytic primitives. |
| `surface` | Standalone EDTSurf, probe 1.4 A, reference radii, voxel atom ownership, distance transform, marching cubes, one smoothing pass, and reference normals. Quality low/medium/high uses detail 3/6/10. |

The EDTSurf wrapper corrects an uninitialized smoothing flag and a radius-index
mismatch that excluded phosphorus. Its original permission notice and local
changes are recorded in [the vendored-source notice](../native/edtsurf/README.md).
CueMol and its source tree are not runtime dependencies. Building a source
distribution requires a C++17 compiler and pybind11; wheels contain the extension.

Ribbon junctions use residue-local parameter tables, analytic scale derivatives,
flat arrow shoulders, fixed section shapes, and reference terminal caps. Remaining
geometric differences include mesh tessellation, chain-break and alternate-location
choices, and modified-base topology. Ribbon junctions share a boundary plane.
Cartoon elements are fitted separately as in
Ribbon2Renderer, with sheet tips kept narrow at direct helix junctions. Maps and
labels remain native; source atoms, colors, bonds, and secondary structure are
unchanged.

## Material and outline audit

The PBR values follow `UmbreonDisplayContext.cpp`. Ambient is evaluated against
unit ambient light. The normalized key direction is `(1,1,1)`, intensity 0.52;
the camera-axis fill is 0.78 with no specular highlight. GGX distribution,
correlated Smith masking, and Schlick Fresnel are shared by GPU/native samples.

| Material | Ambient | Diffuse | Metallic | Roughness | Specular | Reflection |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `default` | .20 | .80 | 0 | .3742032 | .40 | 0 |
| `matte` | .30 | .80 | 0 | .5 | 0 | 0 |
| `diff_metal` | .35 | .30 | 1 | .5491005 | .80 | .10 |
| `spec_metal` | .15 | .60 | 1 | .3742032 | .80 | .65 |
| `metallic_chrome` | .20 | .80 | 1 | .05 | .50 | 0 |
| `metallic_copper` | .20 | .80 | 1 | .15 | .50 | 0 |
| `stone35` | .20 | .80 | 0 | .85 | .25 | 0 |
| `wood31`, `wood14scl2` | .20 | .80 | 0 | .45 | .50 | 0 |

Metal pigments remain the molecular colors. Reflection uses the background;
zero explicit reflection falls back to the PBR F0 environment term. Wood and
stone use the current backend's material intent, without legacy POV-Ray patterns.
`toon1` uses diffuse .8 with brilliance 0; `toon2` uses ambient .3, diffuse .5,
and Phong 10000 with size 50. `shadow` is ambient .75 and `nolighting` is 1.
The toon and metal pairs therefore have distinct appearances. Global
illumination, ambient occlusion, cast shadows, and scene inter-reflections are
not implemented.

`thin`, `normal`, and `thick` widths are .03/.06/.15 A. Lines have a minimum
full width of one output pixel, visibility testing, and depth fog. `outline`
includes self-occlusion contours; `silhouette` suppresses interior contours.
Crease lines are disabled. The screen classifier uses a 12-pixel depth gap,
weak/strong contour hysteresis, normal continuity, and a mesh segment probe to
reject connected folds. Chains shorter than four output pixels are filtered;
continuous junction bars are reconnected before two Chaikin smoothing passes
and .4-output-pixel simplification. Round bands use outside alignment with a
half-output-pixel inner pad. Dedicated GPU PNG, transparent samples, and dedicated ray use
this shared contour path, independent of the triangles' internal diagonals.
Interactive contours use local depth and tangent-plane continuity instead of
tracing and smoothing full image-space chains.

Richardson uses paper `#F0ECDD`, angles 55/-35/80, thresholds .92/.62/.34,
ink scales 1/.74/.38, pressure floor .4, and an absolute display-luma contrast
minimum .15. At 3x, the effective pitch is 2/3 output pixel and full width is
.45 pixel. Stroke length/gap are 50/5 pixels, width jitter .45, length jitter
.5, taper .35, angle jitter 5 degrees, and paper tooth .15 at scale 3.
Hash/noise and stroke coverage have independent numerical reference tests.
Layers are multiplied at each supersample before averaging. The tone recipe
uses diffuse .85, ambient .05, wrap .5, rim power 3.5, rim bias .35, white point
1.2, gamma 2.4, and a highlight knee from .81 to .86. On contexts supporting
`GL_EXT_gpu_shader4`, interactive GLSL evaluates the reference integer hash,
noise, stroke envelopes, and per-sample layer multiplication on the same 3x
grid. Older contexts use approximate floating-point noise while retaining the
paper, angles, thresholds, ink scales, and tone recipe. Live contour joins remain
an approximation of the export's traced and smoothed paths.
Retained standard-ray geometry uses an averaged pencil approximation.

## Validation

```sh
pixi install --locked
pixi run check
uv run --no-project --python .pixi/envs/default/bin/python python -m pytest tests -q
uv run --no-project --python .pixi/envs/default/bin/python python tests/check_cuemol_style.py --gui --benchmark --output .cache/validation
uv build --python .pixi/envs/default/bin/python
```

For image comparison, supply a protein PDB and explicit paths to a compatible
CueMol Node module and its configuration. These are validation-only inputs.

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/audit_appearance.py \
  --structure structure.pdb --output .cache/appearance \
  --reference-module "$CUEMOL_MODULE" --reference-config "$CUEMOL_CONFIG"
uv run --no-project --python .pixi/envs/default/bin/python python tests/compare_appearance.py .cache/appearance
```

The audit covers all 26 profiles and explicit sticks (eight geometries). It
aligns exact coordinates, atom colors, secondary structure, orthographic or
perspective cameras, image dimensions, and renderer properties. Richardson is
compared on the same opaque paper background in both engines; other profiles
use white. The plugin itself never changes the user's background. Re-run with
`--angle`, `--zoom`, `--perspective`, `--transparency`, `--live`, or `--ray`.
`--live` captures the interactive framebuffer instead of the precision export.
Runtime version,
module digest, manifests, logs, images, numerical differences, and contact sheets
are stored in the ignored output directory. Foreground IoU includes shading;
it is not a pure geometric accuracy measure. No image registration is applied.

The image audit used CueMol 2.3.13.523 (`aeacb41`) and the definitions above.
The relevant renderer changes from that build to `3173d8a` add picking names;
the direct exporter, material table, EDTSurf, geometry dimensions, and spline
calculations used here are unchanged. The runtime digest is saved separately
from the source revision so the two are not confused.

Representative measured errors are below. MAE is the mean absolute RGB channel
difference on the union of both foreground masks, in 0..255 units; the blurred
column uses a two-pixel Gaussian. CueMol, GPU PNG, and dedicated ray use 3x
sampling. Protein cases use 1CRN at `medium` quality, and the same camera,
colors, and secondary structure.

| Export and condition | Profile | Foreground IoU | MAE | Blurred MAE |
| --- | --- | ---: | ---: | ---: |
| GPU, opaque, 640x480 | `surface` | 1.0000 | 0.06 | 0.04 |
| GPU, opaque, 640x480 | `cartoon` | 0.9987 | 0.56 | 0.30 |
| GPU, opaque, 640x480 | `default` | 0.9998 | 0.25 | 0.16 |
| GPU, opaque, 640x480 | `richardson` | 0.9967 | 0.72 | 0.18 |
| GPU, perspective, rotated -45 degrees, zoom .8, 640x480 | `richardson` | 0.9893 | 3.28 | 0.51 |
| GPU, transparency .25, 320x240 | `surface` | 0.9986 | 0.37 | 0.22 |
| GPU, transparency .25, 320x240 | `richardson` | 0.9889 | 2.29 | 0.41 |
| GPU, transparency .65, rotated 37 degrees, zoom 1.2, 320x240 | `richardson` | 0.9918 | 1.30 | 0.33 |
| Dedicated ray, opaque, 320x240 | `default` | 0.9963 | 0.72 | 0.47 |
| Dedicated ray, opaque, 320x240 | `surface` | 0.9985 | 0.46 | 0.41 |
| Dedicated ray, opaque, 320x240 | `richardson` | 0.9898 | 2.47 | 0.52 |
| Dedicated ray, transparency .4, 320x240 | `richardson` | 0.9913 | 2.18 | 0.42 |
| Dedicated ray, perspective, rotated -45 degrees, zoom .8, 320x240 | `richardson` | 0.9875 | 4.19 | 0.99 |

Contour placement and renderer-group blending use the shared screen pipeline
described above. The table reports the measured image errors for these cameras
and scenes; it does not establish pixel equality for arbitrary inputs.

The interactive path was also compared directly with CueMol 2.3.15.530
(`af9509e`) using the same aligned input. The live contour approximation has
larger differences than precision exports, especially at occluding joins:

| Live condition | Profile | Foreground IoU | MAE | Blurred MAE |
| --- | --- | ---: | ---: | ---: |
| Opaque, 640x480 | `toon1` | 0.9830 | 6.34 | 1.15 |
| Opaque, 640x480 | `toon2` | 0.9562 | 8.16 | 2.15 |
| Opaque, 640x480 | `richardson` | 0.9719 | 9.18 | 1.33 |
| Opaque, 640x480 | `silhouette` | 0.9870 | 5.03 | 0.99 |
| Perspective, rotated -45 degrees, zoom .8, 640x480 | `toon1` | 0.9768 | 9.05 | 1.65 |
| Perspective, rotated -45 degrees, zoom .8, 640x480 | `richardson` | 0.9654 | 11.50 | 2.00 |

Run the additional GUI regression checks with:

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/check_appearance_gui.py --output .cache/appearance-regression
```

The GUI harness checks 2400x1800 export, native/custom transparency, picking,
state mapping, session reload, and failure restoration. The additional regression
suite covers twelve sheet rotation angles, tiled contour continuity, transfer-buffer
restoration, and overlapping transparent groups. The GPU/ray overlap comparison
measured a foreground MAE of 0.23/255. The benchmark uses 500
residues and 100 synthetic states at 1280x720, recording preparation, warmup,
CPU/GPU storage, peak process RSS, rotation FPS, explicit state-switch FPS,
and actual movie draw rate. Rotation/state benchmarks force completed OpenGL
draws instead of counting Qt repaint requests that may be coalesced. Opaque live
styles bypass the CPU contour/pencil path; dense transparent samples and precision
exports remain more expensive. Reproducible images, builds,
and reports stay ignored; published README
gallery PNGs are versioned and regenerated with `tests/render_gallery.py`.

The additional appearance checks compare tiled and untiled live output in both
camera modes, compile the fallback shader, and compare GPU pencil body pixels
against the independently tested native sampler. The latter measured an average
error of 0.154/255, including RGBA8 attachment quantization.

Run a single profile benchmark independently of the gallery checks:

```sh
uv run --no-project --python .pixi/envs/default/bin/python python tests/check_cuemol_style.py --gui --benchmark-only --benchmark-style richardson --output .cache/benchmark-richardson
```

Use `--benchmark-style toon1` or `ribbon` for the same workload with other
materials. The generated report includes preparation time, peak RSS, mesh/native
CGO/GPU storage, and actual rotation and playback rates. Results depend on the
OpenGL driver and scene coverage; retain machine-specific reports in `.cache/`.
