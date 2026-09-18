"""Interactive materials and surface contours without CPU image transfers."""

from ctypes import c_void_p

import numpy as np


class LivePass:
    tile_size = 510
    scale = 3

    def __init__(self):
        self.framebuffer = 0
        self.textures = []
        self.size = None
        self.bytes = 0

    def clear(self):
        from OpenGL import GL as gl
        from OpenGL.GL.EXT import framebuffer_object as fb

        if self.framebuffer:
            fb.glDeleteFramebuffersEXT(1, [self.framebuffer])
        if self.textures:
            gl.glDeleteTextures(self.textures)
        self.__init__()

    def allocate(self, width, height):
        from OpenGL import GL as gl
        from OpenGL.GL.EXT import framebuffer_object as fb

        if self.size == (width, height):
            return
        limit = int(gl.glGetIntegerv(gl.GL_MAX_TEXTURE_SIZE))
        if max(width, height) > limit or width * height * 24 > 256 * 1024**2:
            raise ValueError(
                "Interactive contour framebuffer exceeds 256 MiB or GL limits"
            )
        self.clear()
        try:
            self.framebuffer = int(fb.glGenFramebuffersEXT(1))
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, self.framebuffer)
            for attachment, internal, layout in (
                (gl.GL_COLOR_ATTACHMENT0, gl.GL_RGBA8, gl.GL_RGBA),
                (gl.GL_COLOR_ATTACHMENT1, gl.GL_RGBA32F, gl.GL_RGBA),
                (
                    gl.GL_DEPTH_ATTACHMENT,
                    gl.GL_DEPTH_COMPONENT24,
                    gl.GL_DEPTH_COMPONENT,
                ),
            ):
                texture = int(gl.glGenTextures(1))
                self.textures.append(texture)
                gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
                for key, value in (
                    (gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST),
                    (gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST),
                    (gl.GL_TEXTURE_WRAP_S, gl.GL_CLAMP_TO_EDGE),
                    (gl.GL_TEXTURE_WRAP_T, gl.GL_CLAMP_TO_EDGE),
                ):
                    gl.glTexParameteri(gl.GL_TEXTURE_2D, key, value)
                gl.glTexImage2D(
                    gl.GL_TEXTURE_2D,
                    0,
                    internal,
                    width,
                    height,
                    0,
                    layout,
                    gl.GL_FLOAT,
                    None,
                )
                fb.glFramebufferTexture2DEXT(
                    gl.GL_FRAMEBUFFER, attachment, gl.GL_TEXTURE_2D, texture, 0
                )
            if (
                fb.glCheckFramebufferStatusEXT(gl.GL_FRAMEBUFFER)
                != gl.GL_FRAMEBUFFER_COMPLETE
            ):
                raise RuntimeError("Interactive contour framebuffer is incomplete")
            self.size = (width, height)
            self.bytes = width * height * 24
        except Exception:
            self.clear()
            raise

    @staticmethod
    def bodies(pool, drawing, projection, viewport, scale=1):
        from OpenGL import GL as gl

        program = pool.body_program(drawing, projection)
        gl.glUniform4f(
            gl.glGetUniformLocation(program, "liveViewport"), *map(float, viewport)
        )
        gl.glUniform1f(gl.glGetUniformLocation(program, "liveScale"), scale)
        for piece in drawing.pieces:
            if piece.mesh.opacity < 0.999999:
                continue
            handles, count, _ = pool.buffer(piece)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, handles[0])
            gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, handles[1])
            gl.glVertexPointer(3, gl.GL_FLOAT, 36, c_void_p(0))
            gl.glNormalPointer(gl.GL_FLOAT, 36, c_void_p(12))
            gl.glColorPointer(3, gl.GL_FLOAT, 36, c_void_p(24))
            gl.glDrawElements(gl.GL_TRIANGLES, count, gl.GL_UNSIGNED_INT, c_void_p(0))

    def draw(self, pool, drawing, modelview, projection, viewport):
        from OpenGL import GL as gl

        if drawing.profile.edges == "none" and drawing.profile.material != "richardson":
            self.bodies(pool, drawing, projection, viewport)
            return
        vx, vy, width, height = map(int, viewport)
        world_per_pixel = 2 / (projection[1, 1] * height)
        if abs(projection[3, 3]) < 0.5:
            world_per_pixel *= max(drawing.fog[0], 1e-6)
        pad = max(2, int(np.ceil(drawing.profile.edge_width / world_per_pixel)) + 1)
        tile = self.tile_size
        tw, th = min(width, tile + 2 * pad), min(height, tile + 2 * pad)
        mode = int(gl.glGetIntegerv(gl.GL_MATRIX_MODE))
        scissor_enabled = gl.glIsEnabled(gl.GL_SCISSOR_TEST)
        scissor_box = gl.glGetIntegerv(gl.GL_SCISSOR_BOX)
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glPushMatrix()
        try:
            for y in range(0, height, tile):
                for x in range(0, width, tile):
                    left = min(max(0, x - pad), width - tw)
                    bottom = min(max(0, y - pad), height - th)
                    crop = np.eye(4)
                    crop[0, 0], crop[1, 1] = width / tw, height / th
                    crop[0, 3] = (width - 2 * left - tw) / tw
                    crop[1, 3] = (height - 2 * bottom - th) / th
                    target = crop @ projection
                    gl.glLoadMatrixd(target.T)
                    self.draw_tile(
                        pool,
                        drawing,
                        target,
                        (vx + left, vy + bottom, tw, th),
                        (left, bottom),
                        (width, height),
                        (vx + x, vy + y, min(tile, width - x), min(tile, height - y)),
                    )
        finally:
            gl.glMatrixMode(gl.GL_PROJECTION)
            gl.glPopMatrix()
            gl.glMatrixMode(mode)
            gl.glViewport(*viewport)
            gl.glScissor(*scissor_box)
            if not scissor_enabled:
                gl.glDisable(gl.GL_SCISSOR_TEST)

    def draw_tile(self, pool, drawing, projection, viewport, offset, total, scissor):
        from OpenGL import GL as gl
        from OpenGL.GL.EXT import framebuffer_object as fb

        previous = int(gl.glGetIntegerv(gl.GL_FRAMEBUFFER_BINDING))
        read_buffer = int(gl.glGetIntegerv(gl.GL_READ_BUFFER))
        draw_buffer = int(gl.glGetIntegerv(gl.GL_DRAW_BUFFER))
        active = int(gl.glGetIntegerv(gl.GL_ACTIVE_TEXTURE))
        unpack = int(gl.glGetIntegerv(gl.GL_PIXEL_UNPACK_BUFFER_BINDING))
        width, height = (int(value) * self.scale for value in viewport[2:])
        try:
            gl.glActiveTexture(gl.GL_TEXTURE0)
            gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, 0)
            self.allocate(width, height)
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, self.framebuffer)
            gl.glDrawBuffers(2, [gl.GL_COLOR_ATTACHMENT0, gl.GL_COLOR_ATTACHMENT1])
            gl.glViewport(0, 0, width, height)
            gl.glDisable(gl.GL_SCISSOR_TEST)
            gl.glClearColor(0, 0, 0, 0)
            gl.glClearDepth(1)
            gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
            self.bodies(
                pool,
                drawing,
                projection,
                (
                    -offset[0] * self.scale,
                    -offset[1] * self.scale,
                    total[0] * self.scale,
                    total[1] * self.scale,
                ),
                self.scale,
            )
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, previous)
            gl.glReadBuffer(read_buffer)
            gl.glDrawBuffer(draw_buffer)
            gl.glViewport(*viewport)
            gl.glEnable(gl.GL_SCISSOR_TEST)
            gl.glScissor(*scissor)
            program = pool.program("live")
            gl.glUseProgram(program)
            for slot, (texture, name) in enumerate(
                zip(self.textures, ("image", "surface", "depth"))
            ):
                gl.glActiveTexture(gl.GL_TEXTURE0 + slot)
                gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
                gl.glUniform1i(gl.glGetUniformLocation(program, name), slot)
            gl.glUniform2f(gl.glGetUniformLocation(program, "imageSize"), width, height)
            gl.glUniform1f(gl.glGetUniformLocation(program, "sampleScale"), self.scale)
            gl.glUniform1i(
                gl.glGetUniformLocation(program, "drawEdges"),
                drawing.profile.edges != "none",
            )
            gl.glUniform2f(
                gl.glGetUniformLocation(program, "viewportOrigin"),
                *map(float, viewport[:2]),
            )
            gl.glUniformMatrix4fv(
                gl.glGetUniformLocation(program, "projection"),
                1,
                True,
                projection.astype(np.float32),
            )
            gl.glUniform1i(
                gl.glGetUniformLocation(program, "outerOnly"),
                drawing.profile.edges == "silhouette",
            )
            gl.glUniform1f(
                gl.glGetUniformLocation(program, "edgeWidth"),
                drawing.profile.edge_width,
            )
            gl.glUniform3f(
                gl.glGetUniformLocation(program, "edgeColor"), *drawing.edge_color
            )
            gl.glUniform3f(
                gl.glGetUniformLocation(program, "background"), *drawing.background
            )
            gl.glUniform2f(gl.glGetUniformLocation(program, "fogRange"), *drawing.fog)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
            gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, 0)
            gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
            gl.glDisableClientState(gl.GL_COLOR_ARRAY)
            gl.glVertexPointer(
                2,
                gl.GL_FLOAT,
                0,
                np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], np.float32),
            )
            gl.glEnable(gl.GL_BLEND)
            gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE_MINUS_SRC_ALPHA)
            gl.glDrawArrays(gl.GL_QUADS, 0, 4)
            gl.glDisable(gl.GL_BLEND)
            gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
            gl.glEnableClientState(gl.GL_COLOR_ARRAY)
        finally:
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, previous)
            gl.glReadBuffer(read_buffer)
            gl.glDrawBuffer(draw_buffer)
            gl.glViewport(*viewport)
            gl.glActiveTexture(active)
            gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, unpack)
