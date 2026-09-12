"""Build the standalone surface and rendering helpers."""

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

setup(
    ext_modules=[
        Pybind11Extension(
            "cuemol_style_in_pymol._edtsurf",
            [
                "native/surface.cpp",
                "native/pencil.cpp",
                "native/visibility.cpp",
                "native/contours.cpp",
                "native/edtsurf/ProteinSurface.cpp",
            ],
            cxx_std=17,
        )
    ],
    cmdclass={"build_ext": build_ext},
)
