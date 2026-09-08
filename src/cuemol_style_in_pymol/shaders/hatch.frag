#version 120
uniform vec2 hatchProjection;
uniform int hatchPerspective;
varying vec3 normalEye;
varying vec3 positionEye;
varying vec3 baseColor;

float randomValue(vec2 key) {
    return fract(sin(dot(key, vec2(12.9898, 78.233))) * 43758.5453);
}

float noise1(float position, float seed) {
    float cell = floor(position);
    float fraction = fract(position);
    fraction *= fraction * (3.0 - 2.0 * fraction);
    return 2.0 * mix(randomValue(vec2(cell, seed)),
        randomValue(vec2(cell + 1.0, seed)), fraction) - 1.0;
}

float noise2(vec2 position, float seed) {
    vec2 cell = floor(position);
    vec2 fraction = fract(position);
    fraction *= fraction * (3.0 - 2.0 * fraction);
    return 2.0 * mix(mix(randomValue(cell + seed),
        randomValue(cell + vec2(1.0, 0.0) + seed), fraction.x),
        mix(randomValue(cell + vec2(0.0, 1.0) + seed),
        randomValue(cell + vec2(1.0) + seed), fraction.x), fraction.y) - 1.0;
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

vec3 bandCoverage(vec3 distance, float radius) {
    // CueMol's default 3x ink sampling, in output-pixel units.
    float halfAA = 0.55 / 3.0;
    return clamp((min(distance + halfAA, vec3(radius))
        - max(distance - halfAA, vec3(-radius))) / (2.0 * halfAA), 0.0, 1.0);
}

float pencil(vec2 uv, vec2 direction, float threshold, float tone, float seed) {
    if (tone >= threshold) return 0.0;
    // The reference's 0.5 px pitch is clamped to 2 device pixels at 3x.
    float pitch = 2.0 / 3.0;
    float along = dot(uv, direction);
    float across = dot(uv, vec2(-direction.y, direction.x));
    float nearest = floor(across / pitch + 0.5);
    float growth = clamp((threshold - tone) * 10.0, 0.0, 1.0);
    vec3 samples = across + direction.y * vec3(-1.0 / 3.0, 0.0, 1.0 / 3.0);
    vec3 coverage0 = vec3(0.0), coverage1 = vec3(0.0), coverage2 = vec3(0.0);
    float drift = noise2(uv / 125.0, seed + 91.0);
    // Include the maximum width, wobble, jitter, stroke tilt, and AA reach.
    for (int offset = -7; offset <= 7; ++offset) {
        float row = nearest + float(offset);
        float phase = randomValue(vec2(row, seed + 1.0));
        float shifted = along + 55.0 * phase;
        float stroke = floor(shifted / 55.0);
        float progress = mod(shifted, 55.0);
        vec2 key = vec2(row + 13.0 * stroke, seed + 7.0);
        float length = min(55.0, 50.0 * (1.0 + 0.5 * (randomValue(key) - 0.5)));
        if (progress >= length) continue;
        float center = row * pitch
            + (randomValue(vec2(row, seed + 3.0)) - 0.5) * 0.44 * pitch;
        center += pitch * noise1(along / 60.0, row + 17.0 * seed);
        float tilt = clamp(1.2 * drift
            + 0.45 * (2.0 * randomValue(key + 31.0) - 1.0), -1.0, 1.0);
        center += 0.08748866 * tilt * (progress - 0.5 * length);
        if (abs(across - center) > 1.1) continue;
        float pressure = 1.0 - 0.45 * randomValue(key + 5.0);
        float radius = 0.225 * growth * pressure
            * (1.0 + 0.45 * noise1(along / 60.0, row + 17.0 * seed + 8.0));
        float entry = 0.35 * (0.35 + 0.35 * randomValue(key + 11.0));
        float release = min(0.60, 0.35 * (0.9 + 1.1 * randomValue(key + 19.0)));
        float fraction = progress / length;
        if (randomValue(key + 23.0) < 0.5) fraction = 1.0 - fraction;
        float head = min(1.0, fraction / entry);
        float tail = min(1.0, (1.0 - fraction) / release);
        radius *= head * (2.0 - head) * tail * (2.0 - tail);
        radius *= 1.0 + 0.27 * noise1(progress / max(8.0, length * 0.4),
            row + 13.0 * stroke + seed + 37.0);
        float markScale = 0.6 + 0.4 * pressure;
        vec3 distance = samples - center;
        coverage0 = max(coverage0, bandCoverage(abs(distance - direction.x / 3.0), radius) * markScale);
        coverage1 = max(coverage1, bandCoverage(abs(distance), radius) * markScale);
        coverage2 = max(coverage2, bandCoverage(abs(distance + direction.x / 3.0), radius) * markScale);
    }
    float tooth = 1.0 - 0.15 * (0.5 + 0.5 * noise2(uv / 3.0, seed + 51.0));
    // Slowly varying stroke envelopes are shared across the nine samples.
    return dot(coverage0 + coverage1 + coverage2, vec3(1.0 / 9.0)) * tooth;
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
