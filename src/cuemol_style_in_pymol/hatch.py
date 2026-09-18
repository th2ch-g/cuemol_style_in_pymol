"""Supersampled pencil compositing in an isolated OpenGL framebuffer."""

from ctypes import c_void_p

import numpy as np

from ._edtsurf import pencil
from .materials import drawing_tone


class HatchPass:
    tile_size = 510

    def __init__(self):
        self.framebuffer = 0
        self.textures = []
        self.size = None
        self.key = None
        self.bytes = 0
        self.contour_key = None
        self.contours = None

    def clear(self):
        from OpenGL import GL as gl
        from OpenGL.GL.EXT.framebuffer_object import glDeleteFramebuffersEXT

        if self.framebuffer:
            glDeleteFramebuffersEXT(1, [self.framebuffer])
        if self.textures:
            gl.glDeleteTextures(self.textures)
        self.__init__()

    def allocate(self, width, height):
        from OpenGL import GL as gl
        from OpenGL.GL.EXT import framebuffer_object as fb

        if self.size == (width, height):
            return
        limit = int(gl.glGetIntegerv(gl.GL_MAX_TEXTURE_SIZE))
        if max(width, height) > limit:
            raise RuntimeError(
                f"The 3x pencil framebuffer exceeds the OpenGL texture limit ({limit})"
            )
        contour_key, contours = self.contour_key, self.contours
        self.clear()
        self.contour_key, self.contours = contour_key, contours
        self.framebuffer = int(fb.glGenFramebuffersEXT(1))
        fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, self.framebuffer)
        self.textures = [int(gl.glGenTextures(1)) for _ in range(5)]
        for index, texture in enumerate(self.textures):
            gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
            gl.glTexParameteri(
                gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST
            )
            gl.glTexParameteri(
                gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST
            )
            gl.glTexParameteri(
                gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_S, gl.GL_CLAMP_TO_EDGE
            )
            gl.glTexParameteri(
                gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_T, gl.GL_CLAMP_TO_EDGE
            )
            depth = index in (2, 4)
            gl.glTexImage2D(
                gl.GL_TEXTURE_2D,
                0,
                gl.GL_DEPTH_COMPONENT24 if depth else gl.GL_RGBA32F,
                width,
                height,
                0,
                gl.GL_DEPTH_COMPONENT if depth else gl.GL_RGBA,
                gl.GL_FLOAT,
                None,
            )
            if index < 3:
                attachment = (
                    gl.GL_DEPTH_ATTACHMENT if depth else gl.GL_COLOR_ATTACHMENT0 + index
                )
                fb.glFramebufferTexture2DEXT(
                    gl.GL_FRAMEBUFFER, attachment, gl.GL_TEXTURE_2D, texture, 0
                )
        if (
            fb.glCheckFramebufferStatusEXT(gl.GL_FRAMEBUFFER)
            != gl.GL_FRAMEBUFFER_COMPLETE
        ):
            raise RuntimeError("The pencil framebuffer is incomplete")
        self.size = (width, height)
        self.bytes = width * height * 56

    def draw(self, pool, drawing, modelview, projection, viewport):
        from OpenGL import GL as gl

        contours = None
        meshes = [p.mesh for p in drawing.pieces if p.mesh.opacity >= 0.999999]
        if meshes and drawing.profile.edges != "none":
            from .raster import buffers

            target = viewport.copy().astype(float)
            target[:2] = 0
            target[2:] *= pool.raster_scale / 3
            key = (
                id(drawing),
                modelview.tobytes(),
                projection.tobytes(),
                target.tobytes(),
                drawing.fog,
            )
            if self.contour_key != key:
                _, _, ink = buffers(
                    drawing,
                    meshes,
                    (modelview, projection, target),
                    drawing.sample_budget,
                )
                self.contours = ink[::-1].copy()
                self.contour_key = key
            contours = self.contours
        mode = int(gl.glGetIntegerv(gl.GL_MATRIX_MODE))
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glPushMatrix()
        try:
            vx, vy, width, height = map(int, viewport)
            if drawing.raster_image is not None and drawing.raster_image.shape[:2] != (
                height * pool.raster_scale,
                width * pool.raster_scale,
            ):
                return
            tile = self.tile_size
            for y in range(0, height, tile):
                for x in range(0, width, tile):
                    left, bottom = max(0, x - 1), max(0, y - 1)
                    right, top = min(width, x + tile + 1), min(height, y + tile + 1)
                    crop = np.eye(4)
                    crop[0, 0], crop[1, 1] = (
                        width / (right - left),
                        height / (top - bottom),
                    )
                    crop[0, 3] = (width - left - right) / (right - left)
                    crop[1, 3] = (height - bottom - top) / (top - bottom)
                    tile_projection = crop @ projection
                    gl.glLoadMatrixd(tile_projection.T)
                    self.draw_tile(
                        pool,
                        drawing,
                        modelview,
                        tile_projection,
                        (vx + left, vy + bottom, right - left, top - bottom),
                        (vx + x, vy + y, min(tile, width - x), min(tile, height - y)),
                        None
                        if contours is None
                        else contours[
                            bottom * pool.raster_scale : top * pool.raster_scale,
                            left * pool.raster_scale : right * pool.raster_scale,
                        ],
                        None
                        if drawing.raster_image is None
                        else (
                            drawing.raster_image[
                                bottom * pool.raster_scale : top * pool.raster_scale,
                                left * pool.raster_scale : right * pool.raster_scale,
                            ],
                            drawing.raster_depth[
                                bottom * pool.raster_scale : top * pool.raster_scale,
                                left * pool.raster_scale : right * pool.raster_scale,
                            ],
                        ),
                    )
            if drawing.raster_image is not None:
                drawing.raster_draws += 1
        finally:
            gl.glMatrixMode(gl.GL_PROJECTION)
            gl.glPopMatrix()
            gl.glMatrixMode(mode)

    def draw_tile(
        self,
        pool,
        drawing,
        modelview,
        projection,
        viewport,
        scissor,
        ink=None,
        samples=None,
    ):
        from OpenGL import GL as gl
        from OpenGL.GL.EXT import framebuffer_object as fb

        previous = int(gl.glGetIntegerv(gl.GL_FRAMEBUFFER_BINDING))
        read_buffer = int(gl.glGetIntegerv(gl.GL_READ_BUFFER))
        draw_buffer = int(gl.glGetIntegerv(gl.GL_DRAW_BUFFER))
        active_texture = int(gl.glGetIntegerv(gl.GL_ACTIVE_TEXTURE))
        pack_buffer = int(gl.glGetIntegerv(gl.GL_PIXEL_PACK_BUFFER_BINDING))
        unpack_buffer = int(gl.glGetIntegerv(gl.GL_PIXEL_UNPACK_BUFFER_BINDING))
        width, height = (
            int(viewport[2]) * pool.raster_scale,
            int(viewport[3]) * pool.raster_scale,
        )
        hatching = drawing.profile.material == "richardson"
        key = (
            id(drawing),
            modelview.tobytes(),
            projection.tobytes(),
            width,
            height,
            drawing.fog,
            drawing.background,
            tuple(a.tobytes() for a in drawing.sample_matrices)
            if drawing.sample_matrices is not None
            else None,
        )
        try:
            gl.glActiveTexture(gl.GL_TEXTURE0)
            gl.glBindBuffer(gl.GL_PIXEL_PACK_BUFFER, 0)
            gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, 0)
            for setting in (
                gl.GL_PACK_ROW_LENGTH,
                gl.GL_PACK_SKIP_ROWS,
                gl.GL_PACK_SKIP_PIXELS,
                gl.GL_UNPACK_ROW_LENGTH,
                gl.GL_UNPACK_SKIP_ROWS,
                gl.GL_UNPACK_SKIP_PIXELS,
            ):
                gl.glPixelStorei(setting, 0)
            self.allocate(width, height)
            if samples is not None and self.key != key:
                gl.glPixelStorei(gl.GL_UNPACK_ALIGNMENT, 1)
                for texture, layout, data in (
                    (3, gl.GL_RGBA, samples[0]),
                    (4, gl.GL_DEPTH_COMPONENT, samples[1]),
                ):
                    gl.glBindTexture(gl.GL_TEXTURE_2D, self.textures[texture])
                    gl.glTexSubImage2D(
                        gl.GL_TEXTURE_2D,
                        0,
                        0,
                        0,
                        width,
                        height,
                        layout,
                        gl.GL_FLOAT,
                        np.ascontiguousarray(data),
                    )
                self.key = key
            elif self.key != key:
                fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, self.framebuffer)
                if hatching:
                    gl.glDrawBuffers(
                        2, [gl.GL_COLOR_ATTACHMENT0, gl.GL_COLOR_ATTACHMENT1]
                    )
                else:
                    gl.glDrawBuffer(gl.GL_COLOR_ATTACHMENT0)
                gl.glViewport(0, 0, width, height)
                gl.glDisable(gl.GL_SCISSOR_TEST)
                gl.glClearColor(0, 0, 0, 0)
                gl.glClearDepth(1)
                gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
                if hatching:
                    gl.glUseProgram(pool.program("hatch"))
                else:
                    pool.body_program(drawing, projection)
                for piece in drawing.pieces:
                    opaque = piece.mesh.opacity >= 0.999999
                    if not opaque:
                        continue
                    gl.glColorMask(opaque, opaque, opaque, opaque)
                    handles, count, _ = pool.buffer(piece)
                    gl.glBindBuffer(gl.GL_ARRAY_BUFFER, handles[0])
                    gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, handles[1])
                    gl.glVertexPointer(3, gl.GL_FLOAT, 36, c_void_p(0))
                    gl.glNormalPointer(gl.GL_FLOAT, 36, c_void_p(12))
                    gl.glColorPointer(3, gl.GL_FLOAT, 36, c_void_p(24))
                    gl.glDrawElements(
                        gl.GL_TRIANGLES, count, gl.GL_UNSIGNED_INT, c_void_p(0)
                    )
                gl.glColorMask(True, True, True, True)
                ink_mask = None if ink is None else ink > 0
                if ink is not None:
                    factor = np.clip(
                        (drawing.fog[1] - ink[ink_mask])
                        / max(drawing.fog[1] - drawing.fog[0], 1e-8),
                        0,
                        1,
                    )
                    background = np.asarray(drawing.background)
                    ink_color = (
                        background
                        + (np.asarray(drawing.edge_color) - background)
                        * factor[:, None]
                    )
                if hatching:
                    gl.glPixelStorei(gl.GL_PACK_ALIGNMENT, 1)
                    gl.glReadBuffer(gl.GL_COLOR_ATTACHMENT0)
                    normals = np.asarray(
                        gl.glReadPixels(0, 0, width, height, gl.GL_RGBA, gl.GL_FLOAT)
                    ).reshape(height, width, 4)
                    gl.glReadBuffer(gl.GL_COLOR_ATTACHMENT1)
                    pigment = np.asarray(
                        gl.glReadPixels(0, 0, width, height, gl.GL_RGBA, gl.GL_FLOAT)
                    ).reshape(height, width, 4)
                    mask = normals[:, :, 3] < 0
                    yy, xx = np.nonzero(mask)
                    uv = np.c_[(xx + 0.5) / 3, (height - yy - 0.5) / 3]
                    if drawing.sample_matrices is not None:
                        target_view, target_projection, target_viewport = (
                            drawing.sample_matrices
                        )
                        transform = (
                            target_projection @ target_view @ np.linalg.inv(modelview)
                        )
                        eye_z = normals[mask, 3]
                        clip_w = projection[3, 2] * eye_z + projection[3, 3]
                        eye_x = (
                            (2 * (xx + 0.5) / width - 1) * clip_w
                            - projection[0, 2] * eye_z
                            - projection[0, 3]
                        ) / projection[0, 0]
                        eye_y = (
                            (2 * (yy + 0.5) / height - 1) * clip_w
                            - projection[1, 2] * eye_z
                            - projection[1, 3]
                        ) / projection[1, 1]
                        target = (
                            np.c_[eye_x, eye_y, eye_z, np.ones(len(xx))] @ transform.T
                        )
                        ndc = target[:, :2] / target[:, 3:4]
                        uv = np.c_[
                            (ndc[:, 0] + 1) * target_viewport[2] / 2,
                            (1 - ndc[:, 1]) * target_viewport[3] / 2,
                        ]
                    direction = None
                    if abs(projection[3, 3]) < 0.5:
                        direction = np.c_[
                            -((2 * (xx + 0.5) / width - 1) + projection[0, 2])
                            / projection[0, 0],
                            -((2 * (yy + 0.5) / height - 1) + projection[1, 2])
                            / projection[1, 1],
                            np.ones(len(xx)),
                        ]
                    fog = np.clip(
                        (drawing.fog[1] + normals[mask, 3])
                        / max(drawing.fog[1] - drawing.fog[0], 1e-8),
                        0,
                        1,
                    )
                    tone = drawing_tone(normals[mask, :3], direction, fog)
                    result = np.zeros_like(pigment)
                    base = np.empty_like(result[:, :, :3])
                    base[:] = np.array([240, 236, 221]) / 255
                    if ink is not None:
                        base[ink_mask] = ink_color
                        result[ink_mask, :3] = ink_color
                        result[ink_mask, 3] = 1
                    result[mask, :3] = pencil(
                        uv, tone, pigment[mask, :3], False, base[mask]
                    )
                    result[mask, 3] = 1
                elif ink is not None:
                    gl.glPixelStorei(gl.GL_PACK_ALIGNMENT, 1)
                    gl.glReadBuffer(gl.GL_COLOR_ATTACHMENT0)
                    result = np.asarray(
                        gl.glReadPixels(0, 0, width, height, gl.GL_RGBA, gl.GL_FLOAT)
                    ).reshape(height, width, 4)
                    result[ink_mask, :3] = ink_color
                    result[ink_mask, 3] = 1
                if ink is not None:
                    depth = np.asarray(
                        gl.glReadPixels(
                            0, 0, width, height, gl.GL_DEPTH_COMPONENT, gl.GL_FLOAT
                        )
                    ).reshape(height, width)
                    z = -ink[ink_mask]
                    ink_depth = 0.5 + 0.5 * (
                        projection[2, 2] * z + projection[2, 3]
                    ) / (projection[3, 2] * z + projection[3, 3])
                    depth[ink_mask] = np.minimum(depth[ink_mask], ink_depth)
                    gl.glBindTexture(gl.GL_TEXTURE_2D, self.textures[4])
                    gl.glPixelStorei(gl.GL_UNPACK_ALIGNMENT, 1)
                    gl.glTexSubImage2D(
                        gl.GL_TEXTURE_2D,
                        0,
                        0,
                        0,
                        width,
                        height,
                        gl.GL_DEPTH_COMPONENT,
                        gl.GL_FLOAT,
                        depth,
                    )
                if hatching or ink is not None:
                    gl.glBindTexture(gl.GL_TEXTURE_2D, self.textures[3])
                    gl.glPixelStorei(gl.GL_UNPACK_ALIGNMENT, 1)
                    gl.glTexSubImage2D(
                        gl.GL_TEXTURE_2D,
                        0,
                        0,
                        0,
                        width,
                        height,
                        gl.GL_RGBA,
                        gl.GL_FLOAT,
                        result,
                    )
                self.key = key
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, previous)
            gl.glReadBuffer(read_buffer)
            gl.glDrawBuffer(draw_buffer)
            gl.glViewport(*viewport)
            gl.glEnable(gl.GL_SCISSOR_TEST)
            gl.glScissor(*scissor)
            program = pool.program("composite")
            gl.glUseProgram(program)
            for unit, texture, name in (
                (
                    0,
                    self.textures[
                        3 if samples is not None or hatching or ink is not None else 0
                    ],
                    "image",
                ),
                (
                    1,
                    self.textures[4 if samples is not None or ink is not None else 2],
                    "depth",
                ),
            ):
                gl.glActiveTexture(gl.GL_TEXTURE0 + unit)
                gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
                gl.glUniform1i(pool.uniform(program, name), unit)
            gl.glUniform2f(pool.uniform(program, "imageSize"), width, height)
            gl.glUniform1i(pool.uniform(program, "samples"), pool.raster_scale)
            gl.glUniform2f(
                pool.uniform(program, "viewportOrigin"),
                float(viewport[0]),
                float(viewport[1]),
            )
            gl.glEnable(gl.GL_BLEND)
            gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE_MINUS_SRC_ALPHA)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
            gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, 0)
            gl.glVertexPointer(
                2,
                gl.GL_FLOAT,
                0,
                np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], np.float32),
            )
            gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
            gl.glDisableClientState(gl.GL_COLOR_ARRAY)
            gl.glDrawArrays(gl.GL_QUADS, 0, 4)
            gl.glDisable(gl.GL_BLEND)
            gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
            gl.glEnableClientState(gl.GL_COLOR_ARRAY)
        finally:
            fb.glBindFramebufferEXT(gl.GL_FRAMEBUFFER, previous)
            gl.glReadBuffer(read_buffer)
            gl.glDrawBuffer(draw_buffer)
            gl.glViewport(*viewport)
            gl.glActiveTexture(active_texture)
            gl.glBindBuffer(gl.GL_PIXEL_PACK_BUFFER, pack_buffer)
            gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, unpack_buffer)
