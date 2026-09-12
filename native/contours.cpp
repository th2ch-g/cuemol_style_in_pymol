#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

namespace py = pybind11;
using Floats = py::array_t<float, py::array::c_style | py::array::forcecast>;
using Doubles = py::array_t<double, py::array::c_style | py::array::forcecast>;
using Faces = py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast>;

namespace {
py::list trace_paths(py::array_t<int, py::array::c_style | py::array::forcecast> begin,
                     py::array_t<int, py::array::c_style | py::array::forcecast> end,
                     py::array_t<bool, py::array::c_style | py::array::forcecast> active) {
    const auto a = begin.unchecked<2>(), b = end.unchecked<2>();
    const auto enabled = active.unchecked<1>();
    if (a.shape(1) != 2 || b.shape(1) != 2 || a.shape(0) != b.shape(0) ||
        a.shape(0) != enabled.shape(0))
        throw std::invalid_argument("Expected matching contour endpoints and mask");
    using Point = std::array<int, 2>;
    struct Chain {
        std::vector<Point> points;
        std::vector<py::ssize_t> edges;
        int first_degree, last_degree;
    };
    std::vector<Chain> chains;
    {
        py::gil_scoped_release release;
        std::map<Point, std::vector<py::ssize_t>> adjacency;
        std::vector<bool> used(a.shape(0), false);
        auto first = [&](py::ssize_t i) { return Point{a(i, 1), a(i, 0)}; };
        auto last = [&](py::ssize_t i) { return Point{b(i, 1), b(i, 0)}; };
        for (py::ssize_t i = 0; i < a.shape(0); ++i)
            if (enabled(i)) {
                adjacency[first(i)].push_back(i);
                adjacency[last(i)].push_back(i);
            }
        auto walk = [&](Point start, py::ssize_t edge) {
            Chain chain;
            chain.points.push_back(start);
            Point current = start;
            chain.first_degree = int(adjacency.at(start).size());
            while (!used[edge]) {
                used[edge] = true;
                chain.edges.push_back(edge);
                current = first(edge) == current ? last(edge) : first(edge);
                chain.points.push_back(current);
                const auto &next = adjacency.at(current);
                if (current == start || next.size() != 2)
                    break;
                edge = next[0] == edge ? next[1] : next[0];
            }
            chain.last_degree = int(adjacency.at(current).size());
            chains.push_back(std::move(chain));
        };
        for (const auto &[point, edges] : adjacency)
            if (edges.size() != 2)
                for (auto edge : edges)
                    if (!used[edge])
                        walk(point, edge);
        for (py::ssize_t i = 0; i < a.shape(0); ++i)
            if (enabled(i) && !used[i])
                walk(first(i), i);
    }
    py::list result;
    for (const auto &chain : chains) {
        py::array_t<double> points({py::ssize_t(chain.points.size()), py::ssize_t(2)});
        py::array_t<py::ssize_t> ids(chain.edges.size());
        auto xy = points.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < xy.shape(0); ++i) {
            xy(i, 0) = chain.points[i][1];
            xy(i, 1) = chain.points[i][0];
        }
        std::copy(chain.edges.begin(), chain.edges.end(), ids.mutable_data());
        result.append(py::make_tuple(points, ids, chain.first_degree, chain.last_degree));
    }
    return result;
}

using Vec = std::array<double, 3>;
Vec subtract(const Vec &a, const Vec &b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec cross(const Vec &a, const Vec &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dot(const Vec &a, const Vec &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

struct Triangle {
    Vec p, a, b, lo, hi;
};
struct Node {
    Vec lo, hi;
    int begin, end, left = -1, right = -1;
};
class SegmentTree {
    std::vector<Triangle> triangles;
    std::vector<Node> nodes;
    int build(int begin, int end) {
        Node node;
        node.begin = begin;
        node.end = end;
        node.lo.fill(std::numeric_limits<double>::infinity());
        node.hi.fill(-std::numeric_limits<double>::infinity());
        for (int i = begin; i < end; ++i)
            for (int k = 0; k < 3; ++k) {
                node.lo[k] = std::min(node.lo[k], triangles[i].lo[k]);
                node.hi[k] = std::max(node.hi[k], triangles[i].hi[k]);
            }
        const int id = int(nodes.size());
        nodes.push_back(node);
        if (end - begin > 12) {
            int axis = 0;
            for (int k = 1; k < 3; ++k)
                if (node.hi[k] - node.lo[k] > node.hi[axis] - node.lo[axis])
                    axis = k;
            const int middle = (begin + end) / 2;
            std::nth_element(triangles.begin() + begin, triangles.begin() + middle,
                             triangles.begin() + end, [axis](const Triangle &a, const Triangle &b) {
                                 return a.lo[axis] + a.hi[axis] < b.lo[axis] + b.hi[axis];
                             });
            const int left = build(begin, middle), right = build(middle, end);
            nodes[id].left = left;
            nodes[id].right = right;
        }
        return id;
    }
    bool query(int index, const Vec &p, const Vec &direction) const {
        const Node &node = nodes[index];
        double low = 0, high = 1;
        for (int k = 0; k < 3; ++k) {
            if (std::abs(direction[k]) < 1e-14) {
                if (p[k] < node.lo[k] || p[k] > node.hi[k])
                    return false;
            } else {
                double a = (node.lo[k] - p[k]) / direction[k],
                       b = (node.hi[k] - p[k]) / direction[k];
                if (a > b)
                    std::swap(a, b);
                low = std::max(low, a);
                high = std::min(high, b);
                if (low > high)
                    return false;
            }
        }
        if (node.left >= 0)
            return query(node.left, p, direction) || query(node.right, p, direction);
        for (int i = node.begin; i < node.end; ++i) {
            const auto &t = triangles[i];
            const Vec h = cross(direction, t.b);
            const double determinant = dot(t.a, h);
            if (std::abs(determinant) < 1e-12)
                continue;
            const Vec s = subtract(p, t.p);
            const double u = dot(s, h) / determinant;
            if (u < 0 || u > 1)
                continue;
            const Vec q = cross(s, t.a);
            const double v = dot(direction, q) / determinant;
            if (v < 0 || u + v > 1)
                continue;
            const double distance = dot(t.b, q) / determinant;
            if (distance > 1e-6 && distance < 1 - 1e-6)
                return true;
        }
        return false;
    }

  public:
    SegmentTree(Doubles vertices, Faces faces) {
        const auto v = vertices.unchecked<2>();
        const auto f = faces.unchecked<2>();
        if (v.shape(1) != 3 || f.shape(1) != 3)
            throw std::invalid_argument("Expected mesh triangles");
        triangles.reserve(f.shape(0));
        for (py::ssize_t i = 0; i < f.shape(0); ++i) {
            Vec a, b, c;
            for (int j = 0; j < 3; ++j)
                if (f(i, j) >= v.shape(0))
                    throw std::invalid_argument("Invalid contour triangle index");
            for (int k = 0; k < 3; ++k) {
                a[k] = v(f(i, 0), k);
                b[k] = v(f(i, 1), k);
                c[k] = v(f(i, 2), k);
            }
            Triangle t;
            t.p = a;
            t.a = subtract(b, a);
            t.b = subtract(c, a);
            for (int k = 0; k < 3; ++k) {
                t.lo[k] = std::min({a[k], b[k], c[k]});
                t.hi[k] = std::max({a[k], b[k], c[k]});
            }
            triangles.push_back(t);
        }
        if (!triangles.empty())
            build(0, int(triangles.size()));
    }
    bool occluded(const Vec &a, const Vec &b) const {
        return !nodes.empty() && query(0, a, subtract(b, a));
    }
};

py::tuple screen_cracks(Floats image, Doubles projection, Doubles vertices, Faces faces,
                        bool outer_only, double far_limit, bool analytic) {
    const auto s = image.unchecked<3>();
    const auto p = projection.unchecked<2>();
    if (s.shape(2) < 4 || p.shape(0) != 4 || p.shape(1) != 4 || s.shape(0) < 1 || s.shape(1) < 1)
        throw std::invalid_argument("Invalid contour image");
    const int h = int(s.shape(0)), w = int(s.shape(1));
    const bool ortho = std::abs(p(3, 3)) > .5;
    const double scale = 2 / (p(1, 1) * h);
    const SegmentTree tree(vertices, faces);
    struct Crack {
        int x, y, axis, owner;
        std::uint8_t flags;
        float z;
    };
    std::vector<Crack> cracks;
    {
        py::gil_scoped_release release;
        auto valid = [&](int x, int y) { return x >= 0 && y >= 0 && x < w && y < h; };
        auto z_at = [&](int x, int y) { return valid(x, y) ? double(s(y, x, 0)) : 0.0; };
        auto clear = [&](int x, int y, int axis) {
            bool near = false, along = false;
            for (int dy = -3; dy <= 3; ++dy)
                for (int dx = -3; dx <= 3; ++dx) {
                    for (int side = 0; side < 2; ++side) {
                        const int xx = x + dx + (axis ? side : 0), yy = y + dy + (axis ? 0 : side);
                        if (valid(xx, yy) && z_at(xx, yy) <= 0)
                            near = true;
                    }
                }
            if (!near)
                return true;
            for (int d = -3; d <= 3; ++d)
                for (int side = 0; side < 2; ++side) {
                    const int xx = x + (axis ? side : d), yy = y + (axis ? d : side);
                    if (valid(xx, yy) && z_at(xx, yy) <= 0)
                        along = true;
                }
            return along;
        };
        auto eye = [&](int x, int y, double z) {
            const double cw = ortho ? 1 : z;
            return Vec{((2 * (x + .5) / w - 1) * cw + p(0, 2) * z - p(0, 3)) / p(0, 0),
                       ((1 - 2 * (y + .5) / h) * cw + p(1, 2) * z - p(1, 3)) / p(1, 1), -z};
        };
        for (int axis : {1, 0})
            for (int y = 0; y < h - (axis ? 0 : 1); ++y)
                for (int x = 0; x < w - (axis ? 1 : 0); ++x) {
                    const int dx = axis ? 1 : 0, dy = axis ? 0 : 1;
                    const int bx = x + dx, by = y + dy;
                    const double a = z_at(x, y), b = z_at(bx, by);
                    if (a <= 0 && b <= 0)
                        continue;
                    const int owner = a > 0 && (b <= 0 || a <= b) ? 0 : 1;
                    const double near = owner ? b : a;
                    if ((a > 0) != (b > 0)) {
                        cracks.push_back(
                            {x, y, axis, owner, std::uint8_t(1 | (owner ? 8 : 0)), float(near)});
                        continue;
                    }
                    if (outer_only && std::max(a, b) <= far_limit)
                        continue;
                    const double px = scale * (ortho ? 1 : std::min(a, b));
                    const double tolerance = std::max(12 * px, 4e-5 * std::min(a, b));
                    const double oa = z_at(x - dx, y - dy), ob = z_at(bx + dx, by + dy);
                    const double sa = oa > 0 ? std::clamp(a - oa, -300 * px, 300 * px) : 0;
                    const double sb = ob > 0 ? std::clamp(b - ob, -300 * px, 300 * px) : 0;
                    const double gap = std::min(std::abs(b - a - sa), std::abs(a - b - sb));
                    if (gap <= .5 * tolerance)
                        continue;
                    const double step = std::abs(b - a);
                    const double inf = std::numeric_limits<double>::infinity();
                    const double left = !valid(x - dx, y - dy) ? 0
                                        : oa > 0               ? std::abs(a - oa)
                                                               : inf;
                    const double right = !valid(bx + dx, by + dy) ? 0
                                         : ob > 0                 ? std::abs(b - ob)
                                                                  : inf;
                    if (!(step > left && step >= right))
                        continue;
                    const bool ridge = sa < -.25 * px && sb < -.25 * px;
                    bool strong = gap > tolerance && (analytic || !ridge);
                    if (strong && analytic)
                        strong = clear(x, y, axis);
                    if (strong && !analytic) {
                        const int nx = owner ? bx : x, ny = owner ? by : y;
                        const int fx = owner ? x : bx, fy = owner ? y : by;
                        const int sign = owner ? 1 : -1;
                        double recession = -1;
                        for (int distance = 1; distance <= 6; ++distance) {
                            const double other =
                                z_at(nx + sign * dx * distance, ny + sign * dy * distance);
                            if (other <= 0)
                                break;
                            recession = std::max(recession, std::abs(other - near) / distance);
                        }
                        strong = recession >= 0 && step > 250 * std::max(recession, px);
                        if (!strong) {
                            Vec na, nb;
                            for (int k = 0; k < 3; ++k) {
                                na[k] = s(y, x, k + 1);
                                nb[k] = s(by, bx, k + 1);
                            }
                            const double norm = std::sqrt(dot(na, na) * dot(nb, nb));
                            strong = norm > 1e-12 && 1 - dot(na, nb) / norm > .3;
                            if (strong) {
                                const double far = std::max(a, b), eps = .5 * px;
                                const double window =
                                    std::min(std::max(.25 * step, 2 * px), step - 2 * eps);
                                strong = window > 0 && !tree.occluded(eye(nx, ny, far - window),
                                                                      eye(fx, fy, far - eps));
                            }
                        }
                    }
                    if (strong || clear(x, y, axis))
                        cracks.push_back({x, y, axis, owner,
                                          std::uint8_t(3 | (owner ? 8 : 0) | (strong ? 32 : 0) |
                                                       (!strong && ridge ? 64 : 0)),
                                          float(near)});
                }
    }
    const py::ssize_t n = cracks.size();
    py::array_t<int> begin({n, py::ssize_t(2)}), end({n, py::ssize_t(2)});
    py::array_t<float> depths(n), owners({n, py::ssize_t(2)});
    py::array_t<std::uint8_t> flags(n);
    auto a = begin.mutable_unchecked<2>(), b = end.mutable_unchecked<2>();
    auto z = depths.mutable_unchecked<1>();
    auto o = owners.mutable_unchecked<2>();
    auto f = flags.mutable_unchecked<1>();
    for (py::ssize_t i = 0; i < n; ++i) {
        const auto &c = cracks[i];
        a(i, 0) = c.x + (c.axis ? 1 : 0);
        a(i, 1) = c.y + (c.axis ? 0 : 1);
        b(i, 0) = a(i, 0) + (c.axis ? 0 : 1);
        b(i, 1) = a(i, 1) + (c.axis ? 1 : 0);
        o(i, 0) = float(c.x + .5 + (c.axis ? c.owner : 0));
        o(i, 1) = float(c.y + .5 + (c.axis ? 0 : c.owner));
        z(i) = c.z;
        f(i) = c.flags;
    }
    return py::make_tuple(begin, end, flags, depths, owners);
}
} // namespace

void bind_contours(py::module_ &module) {
    module.def("screen_cracks", &screen_cracks);
    module.def("trace_paths", &trace_paths);
}
