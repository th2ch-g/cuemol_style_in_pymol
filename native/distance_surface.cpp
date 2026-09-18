// Independent two-pass signed-distance SES with native, bounded workspaces.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
const char *surface_triangles(int flags);

namespace {
using Vec = std::array<double, 3>;
struct Ball { Vec center; double radius; int owner; };
struct Vertex { std::array<float, 3> point, normal; int owner; };
using Face = std::array<std::uint32_t, 3>;
struct Mesh { std::vector<Vertex> vertices; std::vector<Face> faces; };
constexpr int corners[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
    {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};
constexpr int edges[12][2] = {
    {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6},
    {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}};

struct Grid {
    Vec origin;
    std::array<int, 3> size;
    double spacing;
    std::size_t budget;
    std::vector<float> field;
    std::vector<int> owners;
    std::size_t index(const std::array<int, 3> &p) const {
        return (std::size_t(p[0]) * size[1] + p[1]) * size[2] + p[2];
    }
    void fill(const std::vector<Ball> &balls) {
        std::fill(field.begin(), field.end(), 1e6f);
        std::fill(owners.begin(), owners.end(), -1);
        // Sphere-local boxes avoid scanning every sphere in every grid slab.
        // The min reduction retains the reference input order for color ties.
        for (const auto &ball : balls) {
            std::array<int, 3> lo, hi;
            for (int d = 0; d < 3; ++d) {
                const double c = (ball.center[d] - origin[d]) / spacing;
                const double reach = ball.radius / spacing + 2;
                lo[d] = std::max(0, int(std::floor(c - reach)));
                hi[d] = std::min(size[d] - 1, int(std::ceil(c + reach)));
            }
            for (int x = lo[0]; x <= hi[0]; ++x) {
                const double dx = origin[0] + x * spacing - ball.center[0];
                for (int y = lo[1]; y <= hi[1]; ++y) {
                    const double dy = origin[1] + y * spacing - ball.center[1];
                    const double xy = dx * dx + dy * dy;
                    auto i = index({x, y, lo[2]});
                    for (int z = lo[2]; z <= hi[2]; ++z, ++i) {
                        const double dz = origin[2] + z * spacing - ball.center[2];
                        const double squared = xy + dz * dz;
                        const double reach = double(field[i]) + ball.radius;
                        // Overlapping probe balls rarely improve the current
                        // minimum. Reject them before the square root.
                        if (reach >= 0 && squared >= reach * reach) continue;
                        const float value = float(std::sqrt(squared) - ball.radius);
                        if (value < field[i]) {
                            field[i] = value;
                            owners[i] = ball.owner;
                        }
                    }
                }
            }
        }
    }
    Vec gradient(std::array<int, 3> p) const {
        Vec out;
        for (int d = 0; d < 3; ++d) {
            auto a = p, b = p;
            a[d] = std::max(p[d] - 1, 0);
            b[d] = std::min(p[d] + 1, size[d] - 1);
            out[d] = field[index(b)] - field[index(a)];
        }
        return out;
    }
    Mesh contour() const {
        Mesh mesh;
        std::unordered_map<std::size_t, std::uint32_t> cache;
        for (int x = 0; x < size[0] - 1; ++x)
        for (int y = 0; y < size[1] - 1; ++y)
        for (int z = 0; z < size[2] - 1; ++z) {
            std::array<std::array<int, 3>, 8> p;
            int flags = 0;
            for (int c = 0; c < 8; ++c) {
                p[c] = {x + corners[c][0], y + corners[c][1], z + corners[c][2]};
                if (field[index(p[c])] <= 0) flags |= 1 << c;
            }
            const char *triangles = surface_triangles(flags);
            for (int t = 0; t < 16 && static_cast<signed char>(triangles[t]) >= 0; t += 3) {
                Face face;
                for (int j = 0; j < 3; ++j) {
                    const int edge = triangles[t + j];
                    auto a = p[edges[edge][0]], b = p[edges[edge][1]];
                    if (a > b) std::swap(a, b);
                    const int axis = a[0] != b[0] ? 0 : (a[1] != b[1] ? 1 : 2);
                    const auto ia = index(a), ib = index(b), key = ia * 3 + axis;
                    auto found = cache.find(key);
                    if (found != cache.end()) { face[j] = found->second; continue; }
                    // MSVert stores interpolated grid positions as float32.
                    const float f0 = field[ia], f1 = field[ib];
                    const float fraction = f0 == f1 ? 0.5f : -f0 / (f1 - f0);
                    const auto ga = gradient(a), gb = gradient(b);
                    Vertex v;
                    Vec normal;
                    double length = 0;
                    for (int d = 0; d < 3; ++d) {
                        v.point[d] = float(a[d] + (b[d] - a[d]) * double(fraction));
                        normal[d] = ga[d] + (gb[d] - ga[d]) * fraction;
                        length += normal[d] * normal[d];
                    }
                    length = std::sqrt(length);
                    for (int d = 0; d < 3; ++d)
                        v.normal[d] = float(length > 1e-10 ? normal[d] / length : (d == 2));
                    v.owner = owners[f0 <= f1 ? ia : ib];
                    face[j] = std::uint32_t(mesh.vertices.size());
                    cache.emplace(key, face[j]);
                    mesh.vertices.push_back(v);
                }
                mesh.faces.push_back(face);
            }
            // Include hash nodes, vector capacity, both passes, and compaction.
            if (mesh.vertices.size() * 240 + mesh.faces.size() * 48 > budget)
                throw std::runtime_error("Distance surface mesh exceeds workspace budget");
        }
        return mesh;
    }
};

int root(std::vector<int> &parents, int i) {
    while (parents[i] != i) { parents[i] = parents[parents[i]]; i = parents[i]; }
    return i;
}

py::tuple distance_surface(
    py::array_t<float, py::array::c_style | py::array::forcecast> positions,
    py::array_t<int, py::array::c_style | py::array::forcecast> elements,
    int detail, double probe, std::size_t max_bytes)
{
    auto p = positions.unchecked<2>();
    auto e = elements.unchecked<1>();
    if (p.shape(1) != 3 || p.shape(0) != e.shape(0) || !p.shape(0)
        || p.shape(0) > 10000000 || detail < 1 || detail > 100
        || !std::isfinite(probe) || probe < 0 || probe > 10 || max_bytes < 65536)
        throw std::invalid_argument("Invalid distance surface coordinates, detail, probe, or budget");
    constexpr double radii[] = {1.2, 1.7, 1.55, 1.52, 1.8, 1.8, 1.7};
    std::vector<Ball> atoms;
    Vec low = {1e6, 1e6, 1e6}, high = {-1e6, -1e6, -1e6};
    double radius = 0;
    for (py::ssize_t i = 0; i < p.shape(0); ++i) {
        if (e(i) < 0 || e(i) > 6) throw std::invalid_argument("Invalid element radius index");
        Ball ball{{p(i, 0), p(i, 1), p(i, 2)}, radii[e(i)], int(i)};
        radius = std::max(radius, ball.radius);
        for (int d = 0; d < 3; ++d) {
            if (!std::isfinite(ball.center[d]) || std::abs(ball.center[d]) > 1e6)
                throw std::invalid_argument("Invalid distance surface position");
            low[d] = std::min(low[d], ball.center[d]);
            high[d] = std::max(high[d], ball.center[d]);
        }
        atoms.push_back(ball);
    }
    Mesh result;
    {
        py::gil_scoped_release release;
        Grid grid;
        const auto max_cells = std::min(std::size_t(64000000), max_bytes / 32);
        auto dimensions = [&](double spacing) {
            const double pad = radius + 2 * probe + 3 * spacing;
            double count = 1;
            for (int d = 0; d < 3; ++d)
                count *= std::ceil((high[d] - low[d] + 2 * pad) / spacing) + 1;
            return count;
        };
        double spacing = std::max(0.15, 1.43 / (1 + 0.2 * (detail - 1)));
        if (dimensions(spacing) > max_cells) {
            double a = spacing, b = spacing;
            while (dimensions(b) > max_cells) b *= 2;
            for (int i = 0; i < 64; ++i) {
                double mid = (a + b) / 2;
                if (dimensions(mid) > max_cells) a = mid; else b = mid;
            }
            spacing = b;
        }
        grid.spacing = spacing;
        const double pad = radius + 2 * probe + 3 * spacing;
        std::size_t cells = 1;
        for (int d = 0; d < 3; ++d) {
            grid.origin[d] = low[d] - pad;
            grid.size[d] = int(std::ceil((high[d] - low[d] + 2 * pad) / spacing)) + 1;
            cells *= grid.size[d];
        }
        grid.budget = max_bytes - cells * 8;
        grid.field.resize(cells);
        grid.owners.resize(cells);
        auto balls = atoms;
        for (auto &b : balls) b.radius += probe;
        grid.fill(balls);
        auto sas = grid.contour();
        if (probe > 0) {
            balls.clear();
            balls.reserve(sas.vertices.size());
            for (const auto &v : sas.vertices)
                balls.push_back({{grid.origin[0] + v.point[0] * spacing,
                                  grid.origin[1] + v.point[1] * spacing,
                                  grid.origin[2] + v.point[2] * spacing}, probe, v.owner});
            sas = {};
            grid.fill(balls);
            result = grid.contour();
        } else result = std::move(sas);
        for (auto &v : result.vertices)
            for (int d = 0; d < 3; ++d) {
                v.point[d] = float(grid.origin[d] + v.point[d] * spacing);
                if (probe > 0) v.normal[d] = -v.normal[d];
            }
        if (probe > 0) {
            std::vector<int> parents(result.vertices.size()), keep(parents.size(), -1);
            std::iota(parents.begin(), parents.end(), 0);
            for (const auto &f : result.faces) {
                parents[root(parents, f[1])] = root(parents, f[0]);
                parents[root(parents, f[2])] = root(parents, f[0]);
            }
            std::vector<int> remap(parents.size(), -1);
            Mesh compact;
            for (std::size_t i = 0; i < result.vertices.size(); ++i) {
                int r = root(parents, int(i));
                const auto &v = result.vertices[i];
                if (keep[r] < 0) {
                    keep[r] = 0;
                    for (const auto &atom : atoms) {
                        double squared = 0;
                        for (int d = 0; d < 3; ++d) {
                            double delta = v.point[d] - atom.center[d];
                            squared += delta * delta;
                        }
                        if (std::sqrt(squared) - atom.radius < 1.5 * probe) { keep[r] = 1; break; }
                    }
                }
                if (keep[r]) { remap[i] = int(compact.vertices.size()); compact.vertices.push_back(v); }
            }
            for (const auto &f : result.faces)
                if (remap[f[0]] >= 0)
                    compact.faces.push_back({std::uint32_t(remap[f[0]]), std::uint32_t(remap[f[1]]), std::uint32_t(remap[f[2]])});
            result = std::move(compact);
        }
    }
    py::array_t<float> vertices({py::ssize_t(result.vertices.size()), py::ssize_t(3)}), normals(vertices.request().shape);
    py::array_t<std::uint32_t> faces({py::ssize_t(result.faces.size()), py::ssize_t(3)});
    py::array_t<std::int32_t> owners(result.vertices.size());
    auto v = vertices.mutable_unchecked<2>(), n = normals.mutable_unchecked<2>();
    auto f = faces.mutable_unchecked<2>();
    auto o = owners.mutable_unchecked<1>();
    for (std::size_t i = 0; i < result.vertices.size(); ++i) {
        o(i) = result.vertices[i].owner;
        if (o(i) < 0 || o(i) >= p.shape(0)) throw std::runtime_error("Invalid surface owner");
        for (int d = 0; d < 3; ++d) { v(i,d) = result.vertices[i].point[d]; n(i,d) = result.vertices[i].normal[d]; }
    }
    for (std::size_t i = 0; i < result.faces.size(); ++i)
        for (int d = 0; d < 3; ++d) f(i,d) = result.faces[i][d];
    return py::make_tuple(vertices, normals, faces, owners);
}
}

void bind_distance_surface(py::module_ &module) {
    module.def("distance_surface", &distance_surface, py::arg("positions"), py::arg("elements"),
               py::arg("detail") = 6, py::arg("probe") = 1.4,
               py::arg("max_bytes") = std::size_t(2048) * 1024 * 1024);
    module.attr("distance_reference_revision") = "af9509eb381c7b8aa3663475fd07a43f42bf08c4";
}
