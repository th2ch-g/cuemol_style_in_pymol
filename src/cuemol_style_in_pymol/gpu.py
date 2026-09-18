"""Compatibility-profile GLSL renderer; never calls the PyMOL command API."""

from collections import OrderedDict
from dataclasses import dataclass, field
from importlib.resources import files

import numpy as np

from .materials import finish, material_coefficients, material_id
from .presets import PBR_MATERIALS


@dataclass
class Piece:
    mesh: object
    atoms: tuple


@dataclass
class Drawing:
    pieces: list
    profile: object
    edge_color: tuple
    pool: object
    name: str = ""
    state: int = 1
    matrices: object = None
    error: str = ""
    draws: int = 0
    anchors: object = None
    keys: tuple = ()
    selected: object = None
    background: tuple = (0.0, 0.0, 0.0)
    sample_key: object = None
    sample_matrices: object = None
    sampled_bytes: int = 0
    sample_budget: int = 2048 * 1024**2
    raster_image: object = None
    raster_depth: object = None
    raster_draws: int = 0
    fog: tuple = (0.0, 1e10)
    extent: object = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])

    def __post_init__(self):
        from .live import bounds

        for piece in self.pieces:
            # Build acceleration data before any draw or mouse callback.
            piece.mesh.face_bounds
        self.bounds = (
            bounds(self.pieces)
            if self.profile.edges != "none" or self.profile.material == "richardson"
            else np.empty((0, 8, 4), np.float32)
        )
        arrays = [p.mesh.vertices for p in self.pieces if len(p.mesh.vertices)]
        if arrays:
            self.extent = [
                np.min([a.min(axis=0) for a in arrays], axis=0).tolist(),
                np.max([a.max(axis=0) for a in arrays], axis=0).tolist(),
            ]

    def get_extent(self):
        return self.extent

    def __getstate__(self):
        # Session restore tasks rebuild geometry; GL handles are context-local.
        return {"extent": self.extent}

    def __setstate__(self, state):
        self.extent = state["extent"]
        self.pool = None

    def __call__(self):
        if self.pool is None or self.error:
            return
        try:
            self.pool.draw(self)
            self.draws += 1
        except Exception as exc:  # noqa: BLE001 - GUI callback boundary.
            self.error = str(exc)
            print(
                f" cuemol_style: OpenGL drawing failed: {exc}. Use refresh after correcting the OpenGL context."
            )


class Pool:
    """Bounded VBO cache shared by all states in a managed view."""

    def __init__(self, budget_mb=256):
        self.budget = int(budget_mb * 1024**2)
        self.buffers = OrderedDict()
        self.programs = {}
        self.uniforms = {}
        self.context = None
        self.bytes = 0
        self.show_selection = True
        self.raster_scale = 3
        self.precise = False
        self.interacting = False
        from .hatch import HatchPass
        from .live import LivePass

        self.hatch = HatchPass()
        self.live = LivePass()
        self.preview = LivePass()
        self.preview.scale = 1
        self.preview.tile_size = 1020

    def clear(self):
        from OpenGL import GL as gl
        from OpenGL import contextdata

        current = contextdata.getContext()
        if current == self.context:
            for handles, _, _ in self.buffers.values():
                gl.glDeleteBuffers(len(handles), handles)
            for program in self.programs.values():
                gl.glDeleteProgram(program)
            self.hatch.clear()
            self.live.clear()
            self.preview.clear()
        else:
            self.hatch.__init__()
            self.live.__init__()
            self.preview.__init__()
        self.buffers.clear()
        self.programs.clear()
        self.uniforms.clear()
        self.bytes = 0
        self.context = current

    def program(self, name):
        from OpenGL import GL as gl
        from OpenGL.GL.shaders import compileProgram, compileShader

        if name not in self.programs:
            root = files(__package__).joinpath("shaders")
            fragment = root.joinpath(name + ".frag").read_text()
            if "#include pencil" in fragment:
                reference = (
                    b"GL_EXT_gpu_shader4" in gl.glGetString(gl.GL_EXTENSIONS).split()
                )
                if reference:
                    fragment = fragment.replace(
                        "#version 120",
                        "#version 120\n#extension GL_EXT_gpu_shader4 : require",
                    )
                fragment = fragment.replace(
                    "#include pencil",
                    root.joinpath(
                        "pencil.frag" if reference else "pencil_fallback.frag"
                    ).read_text(),
                )
            self.programs[name] = compileProgram(
                compileShader(
                    root.joinpath(name + ".vert").read_text(), gl.GL_VERTEX_SHADER
                ),
                compileShader(fragment, gl.GL_FRAGMENT_SHADER),
                validate=False,
            )
        return self.programs[name]

    def uniform(self, program, name):
        key = program, name
        if key not in self.uniforms:
            from OpenGL import GL as gl

            self.uniforms[key] = int(gl.glGetUniformLocation(program, name))
        return self.uniforms[key]

    def buffer(self, piece):
        from OpenGL import GL as gl

        key = id(piece.mesh)
        if key in self.buffers:
            self.buffers.move_to_end(key)
            return self.buffers[key]
        mesh = piece.mesh
        data = np.c_[mesh.vertices, mesh.normals, mesh.colors].astype(np.float32)
        arrays = [(gl.GL_ARRAY_BUFFER, data), (gl.GL_ELEMENT_ARRAY_BUFFER, mesh.faces)]
        count = mesh.faces.size
        size = sum(a.nbytes for _, a in arrays)
        if size > self.budget:
            raise RuntimeError(
                "A GPU buffer exceeds 256 MiB; lower quality or select fewer atoms"
            )
        while self.buffers and self.bytes + size > self.budget:
            _, (handles, _, old_size) = self.buffers.popitem(last=False)
            gl.glDeleteBuffers(len(handles), handles)
            self.bytes -= old_size
        handles = [int(gl.glGenBuffers(1)) for _ in arrays]
        try:
            for handle, (target, array) in zip(handles, arrays):
                gl.glBindBuffer(target, handle)
                gl.glBufferData(target, array.nbytes, array, gl.GL_STATIC_DRAW)
        except Exception:
            gl.glDeleteBuffers(len(handles), handles)
            raise
        value = handles, count, size
        self.buffers[key] = value
        self.bytes += size
        return value

    def body_program(self, drawing, projection):
        from OpenGL import GL as gl

        program = self.program("body")
        gl.glUseProgram(program)
        physical = PBR_MATERIALS.get(drawing.profile.material)
        gl.glUniform1i(self.uniform(program, "principled"), physical is not None)
        if physical is not None:
            gl.glUniform4f(self.uniform(program, "pbr"), *physical[2:])
        gl.glUniform1i(
            self.uniform(program, "material"),
            material_id(drawing.profile.material),
        )
        gl.glUniform4f(
            self.uniform(program, "materialLighting"),
            *(
                physical[:2] + (0.0, 0.0)
                if physical
                else material_coefficients(drawing.profile.material)
            ),
        )
        gl.glUniform4f(
            self.uniform(program, "materialFinish"),
            *finish(drawing.profile.material)[4:8],
        )
        gl.glUniform1i(
            self.uniform(program, "perspective"),
            int(abs(projection[3, 3]) < 0.5),
        )
        gl.glUniform3f(self.uniform(program, "background"), *drawing.background)
        gl.glUniform2f(self.uniform(program, "fogRange"), *drawing.fog)
        gl.glUniform1i(
            self.uniform(program, "pencilPreview"),
            drawing.profile.material == "richardson",
        )
        return program

    def draw(self, drawing):
        from OpenGL import GL as gl
        from OpenGL import contextdata

        current = contextdata.getContext()
        if current != self.context:
            self.clear()
        # PushAttrib does not include GLSL programs or VBO bindings.
        old_program = int(gl.glGetIntegerv(gl.GL_CURRENT_PROGRAM))
        old_array = int(gl.glGetIntegerv(gl.GL_ARRAY_BUFFER_BINDING))
        old_element = int(gl.glGetIntegerv(gl.GL_ELEMENT_ARRAY_BUFFER_BINDING))
        old_client_texture = int(gl.glGetIntegerv(gl.GL_CLIENT_ACTIVE_TEXTURE))
        modelview = np.asarray(gl.glGetDoublev(gl.GL_MODELVIEW_MATRIX)).T.copy()
        projection = np.asarray(gl.glGetDoublev(gl.GL_PROJECTION_MATRIX)).T.copy()
        viewport = np.asarray(gl.glGetIntegerv(gl.GL_VIEWPORT)).copy()
        drawing.matrices = (modelview, projection, viewport)
        gl.glPushAttrib(gl.GL_ALL_ATTRIB_BITS)
        gl.glPushClientAttrib(gl.GL_CLIENT_ALL_ATTRIB_BITS)
        try:
            gl.glDisable(gl.GL_LIGHTING)
            gl.glDisable(gl.GL_FOG)
            gl.glDisable(gl.GL_BLEND)
            gl.glDisable(gl.GL_CULL_FACE)
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glDepthFunc(gl.GL_LEQUAL)
            gl.glDepthMask(True)
            gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
            gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
            gl.glEnableClientState(gl.GL_COLOR_ARRAY)
            gl.glClientActiveTexture(gl.GL_TEXTURE0)
            gl.glDisableClientState(gl.GL_TEXTURE_COORD_ARRAY)
            gl.glClientActiveTexture(gl.GL_TEXTURE1)
            gl.glDisableClientState(gl.GL_TEXTURE_COORD_ARRAY)
            gl.glClientActiveTexture(gl.GL_TEXTURE0)
            if drawing.raster_image is not None or any(
                p.mesh.opacity >= 0.999999 for p in drawing.pieces
            ):
                renderer = (
                    self.hatch
                    if self.precise or drawing.raster_image is not None
                    else self.preview
                    if self.interacting
                    else self.live
                )
                renderer.draw(self, drawing, modelview, projection, viewport)
            if (
                self.show_selection
                and drawing.selected is not None
                and len(drawing.selected)
            ):
                # Native selection markers cannot depth-test against a callback
                # surface. Show source atom positions as an explicit overlay.
                program = self.program("body")
                gl.glUseProgram(program)
                gl.glUniform1i(
                    self.uniform(program, "material"),
                    -1,
                )
                gl.glUniform4f(
                    self.uniform(program, "materialLighting"),
                    *material_coefficients("nolighting"),
                )
                gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
                gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, 0)
                gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
                gl.glDisableClientState(gl.GL_COLOR_ARRAY)
                for texture in (gl.GL_TEXTURE0, gl.GL_TEXTURE1):
                    gl.glClientActiveTexture(texture)
                    gl.glDisableClientState(gl.GL_TEXTURE_COORD_ARRAY)
                gl.glVertexPointer(3, gl.GL_FLOAT, 0, drawing.selected)
                gl.glColor3f(1.0, 0.1, 0.8)
                gl.glDisable(gl.GL_DEPTH_TEST)
                gl.glDisable(gl.GL_VERTEX_PROGRAM_POINT_SIZE)
                gl.glPointSize(7.0)
                gl.glDrawArrays(gl.GL_POINTS, 0, len(drawing.selected))
        finally:
            gl.glPopClientAttrib()
            gl.glPopAttrib()
            gl.glUseProgram(old_program)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, old_array)
            gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, old_element)
            gl.glClientActiveTexture(old_client_texture)
