#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "edtsurf/CommonPara.h"
#include "edtsurf/ProteinSurface.h"

namespace py = pybind11;
void bind_pencil(py::module_ &module);
void bind_visibility(py::module_ &module);
void bind_contours(py::module_ &module);

static py::tuple surface(
    py::array_t<float, py::array::c_style | py::array::forcecast> positions,
    py::array_t<int, py::array::c_style | py::array::forcecast> elements,
    int detail, double probe, std::size_t max_bytes)
{
    const auto p = positions.unchecked<2>();
    const auto e = elements.unchecked<1>();
    if (p.shape(1) != 3 || p.shape(0) != e.shape(0) || !p.shape(0))
        throw std::invalid_argument("Expected nonempty (N, 3) coordinates and N elements");
    if (p.shape(0) > std::numeric_limits<int>::max() || detail < 1 || detail > 100
        || !std::isfinite(probe) || probe < 0 || probe > 10)
        throw std::invalid_argument("Invalid surface size, detail, or probe radius");
    std::vector<edtsurf::atom> atoms(p.shape(0));
    for (py::ssize_t i = 0; i < p.shape(0); ++i) {
        if (e(i) < 0 || e(i) > 6)
            throw std::invalid_argument("Element radius index must be in 0..6");
        for (int j = 0; j < 3; ++j)
            if (!std::isfinite(p(i, j)) || std::abs(p(i, j)) > 1e6)
                throw std::invalid_argument("Surface coordinates must be finite and within 1e6 angstroms");
        atoms[i] = {1, static_cast<int>(i), static_cast<char>(e(i)),
                    p(i, 0), p(i, 1), p(i, 2), ' ', 0};
    }
    edtsurf::ProteinSurface mesh;
    const double radii[] = {1.2, 1.7, 1.55, 1.52, 1.8, 1.8, 1.7};
    for (int i = 0; i < edtsurf::ProteinSurface::NO_RAD_TYPES; ++i)
        mesh.rasrad[i] = radii[std::min(i, 6)];
    mesh.proberadius = probe;
    mesh.fixsf = 1.0 + (detail - 1) * 0.2;
    {
        py::gil_scoped_release release;
        const int last = static_cast<int>(atoms.size()) - 1;
        mesh.initpara(0, last, atoms.data(), false, true);
        std::size_t voxels = 1;
        for (int length : {mesh.plength, mesh.pwidth, mesh.pheight}) {
            if (length < 1 || length > std::numeric_limits<short>::max()
                || voxels > max_bytes / 128 / static_cast<std::size_t>(length))
                throw std::runtime_error("EDTSurf grid exceeds the surface workspace budget");
            voxels *= static_cast<std::size_t>(length);
        }
        mesh.fillvoxels(0, last, false, atoms.data(), true);
        mesh.buildboundary();
        mesh.fastdistancemap();
        mesh.marchingcube(4);
        mesh.laplaciansmooth(1);
        mesh.computenorm();
    }
    py::array_t<float> vertices({mesh.vertnumber, 3}), normals({mesh.vertnumber, 3});
    py::array_t<std::uint32_t> faces({mesh.facenumber, 3});
    py::array_t<std::int32_t> owners(mesh.vertnumber);
    auto v = vertices.mutable_unchecked<2>();
    auto n = normals.mutable_unchecked<2>();
    auto f = faces.mutable_unchecked<2>();
    auto o = owners.mutable_unchecked<1>();
    for (int i = 0; i < mesh.vertnumber; ++i) {
        const auto &src = mesh.verts[i];
        if (src.atomid < 0 || src.atomid >= static_cast<int>(atoms.size()))
            throw std::runtime_error("EDTSurf returned an invalid atom owner");
        v(i, 0) = src.x / mesh.scalefactor - mesh.ptran.x;
        v(i, 1) = src.y / mesh.scalefactor - mesh.ptran.y;
        v(i, 2) = src.z / mesh.scalefactor - mesh.ptran.z;
        n(i, 0) = src.pn.x; n(i, 1) = src.pn.y; n(i, 2) = src.pn.z;
        o(i) = src.atomid;
    }
    for (int i = 0; i < mesh.facenumber; ++i) {
        f(i, 0) = mesh.faces[i].a;
        f(i, 1) = mesh.faces[i].b;
        f(i, 2) = mesh.faces[i].c;
    }
    return py::make_tuple(vertices, normals, faces, owners);
}

PYBIND11_MODULE(_edtsurf, module)
{
    bind_pencil(module);
    bind_visibility(module);
    bind_contours(module);
    module.def("surface", &surface, py::arg("positions"), py::arg("elements"),
               py::arg("detail") = 6, py::arg("probe") = 1.4,
               py::arg("max_bytes") = std::size_t(2048) * 1024 * 1024);
    module.attr("reference_revision") = "3173d8af62e211dd37b943ee53b3d6a632e6b5d7";
}
