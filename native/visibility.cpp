#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using Points = py::array_t<double, py::array::c_style | py::array::forcecast>;
using Faces = py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast>;

namespace {
struct Triangle {
    double x[3], y[3], z[3], area;
    template <class P, class F> Triangle(const P &p, const F &f, py::ssize_t i) {
        for (int k = 0; k < 3; ++k) {
            const auto j = f(i, k);
            if (j >= p.shape(0))
                throw std::invalid_argument("Invalid triangle index");
            x[k] = p(j, 0);
            y[k] = p(j, 1);
            z[k] = p(j, 2);
            if (!std::isfinite(x[k]) || !std::isfinite(y[k]) || !std::isfinite(z[k]))
                throw std::invalid_argument("Non-finite projected triangle");
        }
        area = (x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]);
    }
    double depth(double px, double py, bool inside) const {
        if (std::abs(area) < 1e-14)
            return std::numeric_limits<double>::infinity();
        const double b = ((px - x[0]) * (y[2] - y[0]) - (py - y[0]) * (x[2] - x[0])) / area;
        const double c = ((x[1] - x[0]) * (py - y[0]) - (y[1] - y[0]) * (px - x[0])) / area;
        if (inside && (b < 0 || c < 0 || b + c > 1))
            return std::numeric_limits<double>::infinity();
        return z[0] + b * (z[1] - z[0]) + c * (z[2] - z[0]);
    }
};

py::tuple visibility(Points points, Faces faces, int width, int height, bool return_depth) {
    const auto p = points.unchecked<2>();
    const auto f = faces.unchecked<2>();
    if (p.shape(1) != 3 || f.shape(1) != 3 || width < 1 || height < 1 ||
        std::size_t(width) * height > 150000000 || f.shape(0) > INT32_MAX)
        throw std::invalid_argument("Invalid visibility image dimensions");
    py::array_t<double> depth({height, width});
    py::array_t<std::int32_t> ids({height, width});
    py::array_t<bool> used(f.shape(0));
    auto d = depth.mutable_unchecked<2>();
    auto winner = ids.mutable_unchecked<2>();
    auto mask = used.mutable_unchecked<1>();
    std::fill_n(depth.mutable_data(), std::size_t(width) * height,
                std::numeric_limits<double>::infinity());
    std::fill_n(ids.mutable_data(), std::size_t(width) * height, -1);
    std::fill_n(used.mutable_data(), f.shape(0), false);
    {
        py::gil_scoped_release release;
        for (py::ssize_t i = 0; i < f.shape(0); ++i) {
            const Triangle t(p, f, i);
            const int x0 = int(
                std::clamp(std::floor(*std::min_element(t.x, t.x + 3)), 0.0, double(width - 1)));
            const int x1 =
                int(std::clamp(std::ceil(*std::max_element(t.x, t.x + 3)), 0.0, double(width - 1)));
            const int y0 = int(
                std::clamp(std::floor(*std::min_element(t.y, t.y + 3)), 0.0, double(height - 1)));
            const int y1 = int(
                std::clamp(std::ceil(*std::max_element(t.y, t.y + 3)), 0.0, double(height - 1)));
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) {
                    const double z = t.depth(x + .5, y + .5, true);
                    if (z < d(y, x)) {
                        d(y, x) = z;
                        winner(y, x) = int(i);
                    }
                }
        }
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x)
                if (winner(y, x) >= 0)
                    mask(winner(y, x)) = true;
    }
    if (return_depth)
        return py::make_tuple(ids, used, depth);
    return py::make_tuple(ids, used);
}

py::tuple sample_image(Points points, Faces faces, Points attributes, int width, int height) {
    const auto p = points.unchecked<2>();
    const auto f = faces.unchecked<2>();
    const auto a = attributes.unchecked<2>();
    if (p.shape(1) != 4 || f.shape(1) != 3 || a.shape(0) != p.shape(0) || a.shape(1) < 1 ||
        a.shape(1) > 16 || width < 1 || height < 1 || std::size_t(width) * height > 150000000)
        throw std::invalid_argument("Invalid material image dimensions");
    const int channels = int(a.shape(1));
    py::array_t<float> values({height, width, channels});
    py::array_t<float> depth({height, width});
    auto out = values.mutable_unchecked<3>();
    auto d = depth.mutable_unchecked<2>();
    std::fill_n(values.mutable_data(), std::size_t(width) * height * channels, 0.0f);
    std::fill_n(depth.mutable_data(), std::size_t(width) * height,
                std::numeric_limits<float>::infinity());
    {
        py::gil_scoped_release release;
        for (py::ssize_t i = 0; i < f.shape(0); ++i) {
            const Triangle t(p, f, i);
            if (std::abs(t.area) < 1e-14)
                continue;
            const int x0 = int(std::clamp(std::ceil(*std::min_element(t.x, t.x + 3) - .5), 0.0,
                                          double(width - 1)));
            const int x1 = int(std::clamp(std::floor(*std::max_element(t.x, t.x + 3) - .5), 0.0,
                                          double(width - 1)));
            const int y0 = int(std::clamp(std::ceil(*std::min_element(t.y, t.y + 3) - .5), 0.0,
                                          double(height - 1)));
            const int y1 = int(std::clamp(std::floor(*std::max_element(t.y, t.y + 3) - .5), 0.0,
                                          double(height - 1)));
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) {
                    const double b = ((x + .5 - t.x[0]) * (t.y[2] - t.y[0]) -
                                      (y + .5 - t.y[0]) * (t.x[2] - t.x[0])) /
                                     t.area;
                    const double c = ((t.x[1] - t.x[0]) * (y + .5 - t.y[0]) -
                                      (t.y[1] - t.y[0]) * (x + .5 - t.x[0])) /
                                     t.area;
                    if (b < -1e-10 || c < -1e-10 || b + c > 1 + 1e-10)
                        continue;
                    const double z = t.z[0] + b * (t.z[1] - t.z[0]) + c * (t.z[2] - t.z[0]);
                    if (z >= d(y, x))
                        continue;
                    double weights[3] = {(1 - b - c) * p(f(i, 0), 3), b * p(f(i, 1), 3),
                                         c * p(f(i, 2), 3)};
                    const double sum = weights[0] + weights[1] + weights[2];
                    if (sum <= 1e-12)
                        continue;
                    d(y, x) = float(z);
                    for (int k = 0; k < channels; ++k)
                        out(y, x, k) =
                            float((weights[0] * a(f(i, 0), k) + weights[1] * a(f(i, 1), k) +
                                   weights[2] * a(f(i, 2), k)) /
                                  sum);
                }
        }
    }
    return py::make_tuple(values, depth);
}

py::array_t<bool>
visible_samples(Points points, Faces faces, Points original, Faces original_faces,
                py::array_t<std::int32_t, py::array::c_style | py::array::forcecast> image) {
    const auto p = points.unchecked<2>(), source = original.unchecked<2>();
    const auto f = faces.unchecked<2>(), source_faces = original_faces.unchecked<2>();
    const auto ids = image.unchecked<2>();
    if (p.shape(1) != 3 || source.shape(1) != 3 || f.shape(1) != 3 || source_faces.shape(1) != 3 ||
        !ids.shape(0) || !ids.shape(1))
        throw std::invalid_argument("Invalid visibility inputs");
    py::array_t<bool> result(f.shape(0));
    auto keep = result.mutable_unchecked<1>();
    py::gil_scoped_release release;
    for (py::ssize_t i = 0; i < f.shape(0); ++i) {
        const Triangle t(p, f, i);
        const double x = (t.x[0] + t.x[1] + t.x[2]) / 3;
        const double y = (t.y[0] + t.y[1] + t.y[2]) / 3;
        const int winner = ids(int(std::clamp(y, 0.0, double(ids.shape(0) - 1))),
                               int(std::clamp(x, 0.0, double(ids.shape(1) - 1))));
        if (winner >= source_faces.shape(0))
            throw std::invalid_argument("Invalid visibility triangle");
        keep(i) =
            winner < 0 || (t.z[0] + t.z[1] + t.z[2]) / 3 <=
                              Triangle(source, source_faces, winner).depth(x, y, false) + 2e-6;
    }
    return result;
}
} // namespace

void bind_visibility(py::module_ &module) {
    module.def("visibility", &visibility, py::arg("points"), py::arg("faces"), py::arg("width"),
               py::arg("height"), py::arg("return_depth") = false);
    module.def("visible_samples", &visible_samples);
    module.def("sample_image", &sample_image);
}
