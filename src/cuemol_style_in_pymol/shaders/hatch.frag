#version 120
uniform vec2 hatchProjection;
uniform int hatchPerspective;
varying vec3 normalEye;
varying vec3 positionEye;
varying vec3 baseColor;

float randomValue(vec2 key) {
    return fract(sin(dot(key, vec2(12.9898, 78.233))) * 43758.5453);
}

float drawingTone(vec3 n, vec3 viewDirection) {
    vec3 key = normalize(vec3(1.0, 1.0, 1.0));
    float ndv = dot(n, viewDirection);
    float flash = ndv > 0.0 ? clamp((ndv + 0.5) / 1.5, 0.0, 1.0) : 0.0;
    float ndl = dot(n, key);
    float side = ndl > 0.0 ? clamp((ndl + 0.5) / 1.5, 0.0, 1.0) : 0.0;
    float t = 0.05 + 0.85 * 1.302 * (0.6 * flash + 0.4 * side);
    float rim = pow(1.0 - max(dot(n, viewDirection), 0.0), 3.5);
    t *= 1.0 - rim * (1.0 - 0.35 * clamp(t, 0.0, 1.0));
    t = pow(clamp(t / 1.2, 0.0, 1.0), 2.4);
    t = t <= 0.0031308 ? 12.92 * t : 1.055 * pow(t, 1.0 / 2.4) - 0.055;
    return mix(t, 1.0, smoothstep(0.81, 0.86, t));
}

float pencil(vec2 uv, vec2 direction, float threshold, float tone, float seed) {
    if (tone >= threshold) return 0.0;
    // Keep the finest live lattice resolvable without multisample buffers.
    float pitch = 2.0;
    float along = dot(uv, direction);
    float across = dot(uv, vec2(-direction.y, direction.x));
    float nearest = floor(across / pitch + 0.5);
    float growth = clamp((threshold - tone) * 10.0, 0.0, 1.0);
    float coverage = 0.0;
    for (int offset = -2; offset <= 2; ++offset) {
        float row = nearest + float(offset);
        float phase = randomValue(vec2(row, seed + 1.0));
        float shifted = along + 55.0 * phase;
        float stroke = floor(shifted / 55.0);
        float progress = mod(shifted, 55.0);
        float individuality = randomValue(vec2(row + 13.0 * stroke, seed + 7.0));
        float length = 50.0 * (0.5 + 0.5 * individuality);
        float endFade = smoothstep(0.0, 0.18 * length, progress)
            * (1.0 - smoothstep(0.72 * length, length, progress));
        float center = row * pitch + (phase - 0.5) * 0.44 * pitch;
        center += 0.6 * sin(along / 9.55 + phase * 6.283185);
        center += (individuality - 0.5) * 0.087 * (progress - 0.5 * length);
        float pressure = (0.65 + 0.7 * individuality)
            * (0.85 + 0.15 * sin(along / 5.0 + phase * 13.0));
        float radius = 0.36 * growth * pressure * endFade;
        float distance = abs(across - center);
        // Box-filter coverage preserves the energy of sub-pixel strokes.
        float ink = max(0.0, min(distance + 0.5, radius)
            - max(distance - 0.5, -radius));
        coverage = max(coverage, ink);
    }
    float tooth = 0.85 + 0.15 * randomValue(floor(uv / 2.0) + seed);
    return clamp(coverage * tooth, 0.0, 1.0);
}

void main() {
    vec3 n = normalize(normalEye);
    vec3 viewDirection = hatchPerspective == 0 ? vec3(0.0, 0.0, 1.0)
        : normalize(-positionEye);
    if (dot(n, viewDirection) < 0.0) n = -n;
    float tone = drawingTone(n, viewDirection);
    vec2 uv = positionEye.xy * hatchProjection;
    if (hatchPerspective != 0) uv /= max(-positionEye.z, 0.001);
    uv.y = -uv.y;
    vec3 paper = vec3(0.941, 0.925, 0.867);
    vec3 ink = baseColor * (0.4 + 0.6 * tone);
    vec3 luma = vec3(0.2126, 0.7152, 0.0722);
    float target = dot(paper, luma) * 0.85;
    ink *= min(1.0, target / max(dot(ink, luma), 0.0001));
    vec3 result = paper;
    result *= mix(vec3(1.0), ink, pencil(uv, vec2(0.573576, 0.819152), 0.92, tone, 0.0));
    result *= mix(vec3(1.0), ink * 0.74, pencil(uv, vec2(0.819152, -0.573576), 0.62, tone, 1.0));
    result *= mix(vec3(1.0), ink * 0.38, pencil(uv, vec2(0.173648, 0.984808), 0.34, tone, 2.0));
    gl_FragColor = vec4(clamp(result, 0.0, 1.0), 1.0);
}
