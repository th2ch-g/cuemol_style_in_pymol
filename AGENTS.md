# AGENTS.md

`CLAUDE.md` imports this file with `@AGENTS.md`.

## Overview

This repository provides the standalone `cuemol_style` command for PyMOL 3.1.
It implements CueMol-inspired molecular geometry, live GPU materials, picking,
MD state playback, transparency, and PNG/ray export. CueMol and mdtbx are not
runtime dependencies. Read `README.md` and `docs/pymol_cuemol.md` before changing
the public interface or rendering behavior.

## Development

The pixi environment supplies Python 3.10 and PyMOL. Python package dependencies
are declared in `pyproject.toml`; keep `pixi.lock` synchronized. Run Python
through `uv` with the pixi interpreter.

```sh
pixi install --locked
pixi run check
pixi run format
uv run --no-project --python .pixi/envs/default/bin/python python -m pytest tests -q
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --output .cache/headless
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/check_cuemol_style.py --gui --benchmark --output .cache/gui
uv build
```

GUI validation requires Qt and a compatibility OpenGL 2.1 / GLSL 1.20 context.
Use a separate PyMOL process for validation; preserve existing user sessions.
For rendering changes, inspect generated GPU and ray images in addition to
running tests. The benchmark uses 500 residues and 100 synthetic states at
1280 x 720; after preparation, target at least 30 FPS rotation and 15 states/s
movie playback. Report preparation time and memory as well as frame rates.

## Architecture

- `src/cuemol_style_in_pymol/__init__.py`: command registration and dispatch.
- `controller.py`: managed views, refresh/reset, visibility and state lifecycle.
- `source.py`, `geometry.py`, `mesh.py`: source snapshots and molecular meshes.
- `presets.py`, `materials.py`: profiles, palettes, and native material samples.
- `gpu.py`, `shaders/`: isolated OpenGL callbacks, buffers, and live shaders.
- `picking.py`, `export.py`: original-atom selection and image export.
- `tests/test_cuemol_*.py`: geometry and real headless PyMOL checks.
- `tests/check_cuemol_style.py`, `tests/cuemol_benchmark.py`: GUI integration
  and performance validation.

## Rendering constraints

- Do not patch PyMOL source or replace its standard `cmd` functions. Package
  import must not apply styles or change PyMOL settings.
- Preserve source coordinates, bonds, colors, secondary structure, background,
  and unrelated representations. Restore owned visibility and temporary state
  on reset, export failure, deletion, and session reload.
- Prepare source data and all loaded states outside drawing callbacks. During
  a callback, only the existing read-only state/frame queries may enter PyMOL;
  do not add reentrant API calls or scene mutations.
- Restore OpenGL state and the active Qt framebuffer after every callback.
- Opaque bodies use GLSL. Transparent bodies and retained standard-ray geometry
  use native CGO; they cannot evaluate the per-pixel Richardson pencil shader.
  Keep those rendering differences explicit in the guides.
- Compare geometry and material defaults with the reference definitions recorded
  in the guide. Distinguish matched numeric settings from approximate geometry,
  procedural textures, lighting, and image-space effects.

## Documentation and repository conventions

- Reply in concise Japanese. Write code comments and commit messages in English.
- Keep documentation pages in Markdown. Maintain the English guide and its
  Japanese counterpart under `docs/ja/` together.
- Keep paths portable and avoid machine-specific configuration.
- Scope text/file searches to relevant paths and names.
- Ignore reproducible caches, downloaded coordinates, environments, builds,
  and validation output. Published README gallery PNGs are the explicit
  exception: version them so GitHub can display the gallery.
- The gallery belongs between Install and Use. Cover every named profile and
  supported geometry; keep image dimensions, display widths, and comparison
  cameras consistent. Regenerate with `tests/render_gallery.py`.
- Keep changes focused. Append a blank line and
  `Co-authored-by: Codex <noreply@openai.com>` to Git commit messages.
