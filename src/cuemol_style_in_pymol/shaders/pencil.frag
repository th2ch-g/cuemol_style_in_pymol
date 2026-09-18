// Integer hash and stroke geometry matching the native reference sampler.
unsigned int hashMix(unsigned int value) {
    value = (value ^ (value >> 16u)) * 0x7feb352du;
    value = (value ^ (value >> 15u)) * 0x846ca68bu;
    return value ^ (value >> 16u);
}
unsigned int hashCombine(unsigned int left, unsigned int right) {
    return hashMix(left ^ (right + 0x9e3779b9u + (left << 6u) + (left >> 2u)));
}
unsigned int pencilKey(int layer, int row, int column, int stream) {
    unsigned int h = hashMix(unsigned int(layer) ^ 0x9e3779b9u);
    h = hashCombine(h, unsigned int(layer));
    h = hashCombine(h, unsigned int(row));
    h = hashCombine(h, unsigned int(column));
    return hashCombine(h, unsigned int(stream));
}
float random01(unsigned int h) { return float(h >> 8u) / 16777216.0; }
float randomSigned(unsigned int h) { return 2.0 * random01(h) - 1.0; }
float pencilNoise(unsigned int seed, float x) {
    int cell = int(floor(x));
    float f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    return mix(randomSigned(hashMix(hashCombine(seed, unsigned int(cell)))),
               randomSigned(hashMix(hashCombine(seed, unsigned int(cell + 1)))), f);
}
float noiseAt(unsigned int seed, int x, int y) {
    return randomSigned(hashMix(hashCombine(hashCombine(seed, unsigned int(x)), unsigned int(y))));
}
float pencilNoise2(unsigned int seed, vec2 xy) {
    ivec2 cell = ivec2(floor(xy));
    vec2 f = fract(xy);
    f = f * f * (3.0 - 2.0 * f);
    return mix(mix(noiseAt(seed, cell.x, cell.y), noiseAt(seed, cell.x + 1, cell.y), f.x),
               mix(noiseAt(seed, cell.x, cell.y + 1), noiseAt(seed, cell.x + 1, cell.y + 1), f.x), f.y);
}
float pencilLayer(vec2 xy, vec2 direction, float tone, float threshold, float layerId) {
    if (tone >= threshold) return 0.0;
    int layer = int(layerId);
    float across = dot(xy, vec2(-direction.y, direction.x)), along = dot(xy, direction);
    float growth = clamp((threshold - tone) * 10.0, 0.0, 1.0);
    float reach = 0.675 * 1.45 * 1.27 + 0.55 + 0.44 + 2.0 + 0.0874886635 * 165.0 / 2.0;
    int first = int(floor((across - reach) / 2.0));
    int last = int(floor((across + reach) / 2.0));
    float result = 0.0;
    for (int index = 0; index < 13; ++index) {
        int row = first + index;
        if (row > last) break;
        float center = 2.0 * float(row) + (random01(pencilKey(layer, 0, row, 1)) - 0.5) * 0.88;
        center += 2.0 * pencilNoise(pencilKey(layer, 0, row, 3), along / 180.0);
        float radius = 0.675 * growth * (1.0 + 0.45 * pencilNoise(pencilKey(layer, 0, row, 4), along / 180.0));
        float phased = along + 165.0 * random01(pencilKey(layer, 0, row, 5));
        int stroke = int(floor(phased / 165.0));
        float phase = phased - float(stroke) * 165.0;
        float extent = clamp(150.0 * (1.0 + 0.5 * (random01(pencilKey(layer, row, stroke, 9)) - 0.5)), 2.0, 165.0);
        if (phase >= extent) continue;
        float pressure = 1.0 - 0.45 * random01(pencilKey(layer, row, stroke, 10));
        radius *= pressure;
        float strokeCenter = along + (float(stroke) + 0.5) * 165.0 - phased;
        float drift = pencilNoise2(pencilKey(layer, 0, 0, 12), vec2(float(row) * 2.0, strokeCenter) / 375.0);
        float angle = clamp(1.2 * drift + 0.45 * randomSigned(pencilKey(layer, row, stroke, 11)), -1.0, 1.0);
        center += 0.0874886635 * angle * (phase - extent / 2.0);
        unsigned int taper = pencilKey(layer, row, stroke, 13);
        float entry = clamp(0.35 * (0.35 + 0.35 * random01(taper)), 0.02, 0.45);
        float release = clamp(0.35 * (0.9 + 1.1 * random01(hashMix(taper))), 0.05, 0.60);
        if ((taper & 1u) != 0u) { float swap = entry; entry = release; release = swap; }
        float head = clamp(phase / extent / entry, 0.0, 1.0);
        float tail = clamp((1.0 - phase / extent) / release, 0.0, 1.0);
        radius *= head * (2.0 - head) * tail * (2.0 - tail);
        radius *= 1.0 + 0.27 * pencilNoise(pencilKey(layer, row, stroke, 14), phase / max(8.0, extent * 0.4));
        float distance = abs(across - center);
        float fraction = clamp((min(distance + 0.55, radius) - max(distance - 0.55, -radius)) / 1.1, 0.0, 1.0);
        result = max(result, fraction * (0.6 + 0.4 * pressure));
    }
    return result * (1.0 - 0.15 * (0.5 + 0.5 * pencilNoise2(pencilKey(layer, 0, 0, 6), xy / 9.0)));
}
