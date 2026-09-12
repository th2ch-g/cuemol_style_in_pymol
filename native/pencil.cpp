#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace {
using Word = std::uint32_t;
float clamp(float x) { return std::clamp(x, 0.0f, 1.0f); }
Word mix(Word value) {
    value = (value ^ (value >> 16)) * 0x7feb352dU;
    value = (value ^ (value >> 15)) * 0x846ca68bU;
    return value ^ (value >> 16);
}
Word combine(Word left, Word right) {
    return mix(left ^ (right + 0x9e3779b9U + (left << 6) + (left >> 2)));
}
Word key(int layer, int row, int column, int stream) {
    Word hash = mix(Word(layer) ^ 0x9e3779b9U);
    for (int value : {layer, row, column, stream}) hash = combine(hash, Word(value));
    return hash;
}
float random01(Word hash) { return float(hash >> 8) / 16777216.0f; }
float random_signed(Word hash) { return 2 * random01(hash) - 1; }
float lerp(float a, float b, float t) { return a + t * (b - a); }
float ease(float t) { return t * t * (3 - 2 * t); }
float noise(Word seed, float x) {
    const int cell = int(std::floor(x));
    return lerp(random_signed(mix(combine(seed, Word(cell)))),
                random_signed(mix(combine(seed, Word(cell + 1)))), ease(x - cell));
}
float noise(Word seed, float x, float y) {
    const int a = int(std::floor(x)), b = int(std::floor(y));
    auto sample = [seed](int i, int j) {
        return random_signed(mix(combine(combine(seed, Word(i)), Word(j))));
    };
    const float tx = ease(x - a), ty = ease(y - b);
    return lerp(lerp(sample(a, b), sample(a + 1, b), tx),
                lerp(sample(a, b + 1), sample(a + 1, b + 1), tx), ty);
}

struct Pencil {
    int layer;
    float cosine, sine, threshold, darkness;
    float coverage(float x, float y, float tone) const {
        if (tone >= threshold) return 0;
        // Geometry is evaluated on CueMol's default three-sample grid.
        constexpr float pitch = 2, half_filter = 0.55f, length = 150, period = 165;
        const float across = -sine * x + cosine * y;
        const float along = cosine * x + sine * y;
        const float growth = clamp((threshold - tone) * 10);
        constexpr float tilt = 0.0874886635f;
        constexpr float reach = 0.675f * 1.45f * 1.27f + half_filter + 0.44f + 2 + tilt * period / 2;
        const int first = int(std::floor((across - reach) / pitch));
        const int last = int(std::floor((across + reach) / pitch));
        float result = 0;
        for (int row = first; row <= last; ++row) {
            float center = pitch * row + (random01(key(layer, 0, row, 1)) - 0.5f) * 0.88f;
            center += 2 * noise(key(layer, 0, row, 3), along / 180);
            float radius = 0.675f * growth * (1 + 0.45f * noise(key(layer, 0, row, 4), along / 180));
            const float phased = along + period * random01(key(layer, 0, row, 5));
            const int stroke = int(std::floor(phased / period));
            const float phase = phased - stroke * period;
            const float extent = std::clamp(length * (1 + 0.5f * (random01(key(layer, row, stroke, 9)) - 0.5f)), 2.0f, period);
            if (phase >= extent) continue;
            const float pressure = 1 - 0.45f * random01(key(layer, row, stroke, 10));
            radius *= pressure;
            const float stroke_center = along + (stroke + 0.5f) * period - phased;
            const float drift = noise(key(layer, 0, 0, 12), row * pitch / 375, stroke_center / 375);
            const float angle = std::clamp(1.2f * drift + 0.45f * random_signed(key(layer, row, stroke, 11)), -1.0f, 1.0f);
            center += tilt * angle * (phase - extent / 2);
            const Word taper_seed = key(layer, row, stroke, 13);
            float entry = std::clamp(0.35f * (0.35f + 0.35f * random01(taper_seed)), 0.02f, 0.45f);
            float release = std::clamp(0.35f * (0.9f + 1.1f * random01(mix(taper_seed))), 0.05f, 0.60f);
            if (taper_seed & 1) std::swap(entry, release);
            const float head = clamp(phase / extent / entry);
            const float tail = clamp((1 - phase / extent) / release);
            radius *= head * (2 - head) * tail * (2 - tail);
            radius *= 1 + 0.27f * noise(key(layer, row, stroke, 14), phase / std::max(8.0f, extent * 0.4f));
            const float distance = std::abs(across - center);
            const float fraction = clamp((std::min(distance + half_filter, radius) - std::max(distance - half_filter, -radius)) / (2 * half_filter));
            result = std::max(result, fraction * (0.6f + 0.4f * pressure));
        }
        return result * (1 - 0.15f * (0.5f + 0.5f * noise(key(layer, 0, 0, 6), x / 9, y / 9)));
    }
};

const std::array<Pencil, 3> pencils = [] {
    std::array<Pencil, 3> result;
    const float angles[] = {55, -35, 80}, thresholds[] = {0.92f, 0.62f, 0.34f};
    const float darkness[] = {1, 0.74f, 0.38f};
    for (int i = 0; i < 3; ++i) {
        const float angle = angles[i] * 0.017453292519943295f;
        result[i] = {i, std::cos(angle), std::sin(angle), thresholds[i], darkness[i]};
    }
    return result;
}();

std::array<float, 3> pencil_color(float x, float y, float tone, const std::array<float, 3>& albedo,
                          const std::array<float, 3>& paper) {
    const std::array<float, 3> luma = {0.2126f, 0.7152f, 0.0722f};
    std::array<float, 3> ink, result = paper;
    float ink_luma = 0, paper_luma = 0;
    for (int channel = 0; channel < 3; ++channel) {
        ink[channel] = albedo[channel] * (0.4f + 0.6f * tone);
        ink_luma += ink[channel] * luma[channel];
        paper_luma += paper[channel] * luma[channel];
    }
    if (paper_luma - ink_luma < .15f) {
        const float target = paper_luma >= .15f ? paper_luma - .15f : std::min(1.f, paper_luma + .15f);
        for (float &c : ink) c = ink_luma > 1e-4f ? clamp(c * target / ink_luma) : target;
    }
    for (const Pencil &pencil : pencils) {
        const float coverage = pencil.coverage(x, y, tone);
        for (int channel = 0; channel < 3; ++channel)
            result[channel] *= 1 - coverage * (1 - ink[channel] * pencil.darkness);
    }
    return result;
}

py::array_t<float> pencil_samples(
    py::array_t<float, py::array::c_style | py::array::forcecast> coordinates,
    py::array_t<float, py::array::c_style | py::array::forcecast> tones,
    py::array_t<float, py::array::c_style | py::array::forcecast> albedo,
    bool antialias, py::object base)
{
    const auto xy = coordinates.unchecked<2>();
    const auto tone = tones.unchecked<1>();
    const auto colors = albedo.unchecked<2>();
    if (xy.shape(1) != 2 || colors.shape(1) != 3 || xy.shape(0) != tone.shape(0) || xy.shape(0) != colors.shape(0))
        throw std::invalid_argument("Expected coordinates (N,2), tones (N), albedo (N,3)");
    py::array_t<float, py::array::c_style | py::array::forcecast> canvas;
    const float *pixels = nullptr;
    if (!base.is_none()) {
        canvas = decltype(canvas)::ensure(base);
        if (!canvas || canvas.ndim() != 2 || canvas.shape(0) != xy.shape(0) || canvas.shape(1) != 3)
            throw std::invalid_argument("Expected optional base colors (N,3)");
        pixels = canvas.data();
    }
    py::array_t<float> output({xy.shape(0), py::ssize_t(3)});
    auto result = output.mutable_unchecked<2>();
    py::gil_scoped_release release;
    for (py::ssize_t i = 0; i < xy.shape(0); ++i) {
        if (!std::isfinite(xy(i, 0)) || !std::isfinite(xy(i, 1))
            || std::abs(xy(i, 0)) > 1e7 || std::abs(xy(i, 1)) > 1e7
            || !std::isfinite(tone(i)))
            throw std::invalid_argument("Pencil coordinates and tone must be finite and bounded");
        std::array<float, 3> color = {colors(i, 0), colors(i, 1), colors(i, 2)};
        for (float c : color)
            if (!std::isfinite(c)) throw std::invalid_argument("Pencil albedo must be finite");
        std::array<float, 3> paper = {240.f / 255, 236.f / 255, 221.f / 255};
        if (pixels) for (int c = 0; c < 3; ++c) {
            paper[c] = pixels[i * 3 + c];
            if (!std::isfinite(paper[c])) throw std::invalid_argument("Pencil base must be finite");
        }
        std::array<float, 3> sum = {};
        const int radius = antialias ? 1 : 0;
        for (int dy = -radius; dy <= radius; ++dy)
            for (int dx = -radius; dx <= radius; ++dx) {
                auto value = pencil_color(3 * xy(i, 0) + dx, 3 * xy(i, 1) + dy, tone(i), color, paper);
                for (int c = 0; c < 3; ++c) sum[c] += value[c];
            }
        for (int c = 0; c < 3; ++c) result(i, c) = sum[c] / (antialias ? 9 : 1);
    }
    return output;
}
}

void bind_pencil(py::module_ &module) {
    module.def("pencil", &pencil_samples, py::arg("coordinates"), py::arg("tones"),
               py::arg("albedo"), py::arg("antialias") = true, py::arg("base") = py::none());
}
