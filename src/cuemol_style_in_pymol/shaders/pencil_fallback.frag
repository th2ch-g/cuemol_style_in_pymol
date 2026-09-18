float grain(vec2 p) {
    return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453);
}
float noise(vec2 p) {
    vec2 cell = floor(p), f = fract(p);
    f = f * f * (3.0 - 2.0 * f);
    return mix(mix(grain(cell), grain(cell + vec2(1, 0)), f.x),
               mix(grain(cell + vec2(0, 1)), grain(cell + vec2(1, 1)), f.x), f.y);
}
float pencilLayer(vec2 xy, vec2 direction, float tone, float threshold, float seed) {
    seed = seed * 28.0 + 1.0;
    if (tone >= threshold) return 0.0;
    float along = dot(xy, direction), across = dot(xy, vec2(-direction.y, direction.x));
    float growth = clamp((threshold - tone) * 10.0, 0.0, 1.0);
    float filterWidth = max(0.55, fwidth(across) * 0.5);
    float coverage = 0.0;
    float first = floor(across / 2.0);
    for (int i = -2; i <= 2; ++i) {
        float row = first + float(i);
        float center = row * 2.0 + (grain(vec2(row, seed)) - 0.5) * 0.88;
        center += 4.0 * (noise(vec2(along / 180.0, row + seed)) - 0.5);
        float phased = along + 165.0 * grain(vec2(row, seed + 5.0));
        float stroke = floor(phased / 165.0), phase = mod(phased, 165.0);
        float extent = min(165.0, 150.0 * (0.75 + 0.5 * grain(vec2(row + seed, stroke))));
        float pressure = 1.0 - 0.45 * grain(vec2(row + seed + 10.0, stroke));
        float radius = 0.675 * growth * pressure;
        radius *= 0.55 + 0.9 * noise(vec2(along / 180.0, row + seed + 4.0));
        float head = clamp(phase / extent / 0.18, 0.0, 1.0);
        float tail = clamp((1.0 - phase / extent) / 0.45, 0.0, 1.0);
        radius *= head * (2.0 - head) * tail * (2.0 - tail);
        float distance = abs(across - center);
        coverage += clamp((min(distance + filterWidth, radius) - max(distance - filterWidth, -radius))
                          / (2.0 * filterWidth), 0.0, 1.0) * (0.6 + 0.4 * pressure);
    }
    return min(1.0, coverage) * (0.85 + 0.15 * noise(xy / 9.0 + seed));
}
