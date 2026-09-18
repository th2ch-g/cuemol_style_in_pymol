#version 120
uniform sampler2D image;
uniform sampler2D surface;
uniform sampler2D depth;
uniform vec2 imageSize;
uniform float sampleScale;
uniform int drawEdges;
uniform vec2 viewportOrigin;
uniform mat4 projection;
uniform int outerOnly;
uniform float edgeWidth;
uniform vec3 edgeColor;
uniform vec3 background;
uniform vec2 fogRange;
vec3 eyePosition(vec2 uv, float z) {
    float w = projection[2][3] * z + projection[3][3];
    return vec3(((uv * 2.0 - 1.0) * w - vec2(projection[2][0], projection[2][1]) * z
                 - vec2(projection[3][0], projection[3][1]))
                / vec2(projection[0][0], projection[1][1]), z);
}
vec4 sampleColor(vec2 uv, out float resultDepth) {
    vec4 c = texture2D(image, uv);
    vec4 s = texture2D(surface, uv);
    float z = texture2D(depth, uv).r;
    float cameraDepth = abs(projection[3][3]) < 0.5 ? max(fogRange.x, 0.001) : 1.0;
    float pixel = 2.0 * cameraDepth / (projection[1][1] * imageSize.y);
    float radius = max(0.5 * sampleScale, edgeWidth / pixel * 0.5);
    float ink = 0.0;
    float inkZ = s.a;
    vec3 p = eyePosition(uv, s.a);
    if (drawEdges != 0) {
        for (int i = 0; i < 8; ++i) {
            float angle = float(i) * 0.7853981634;
            vec2 offset = vec2(cos(angle), sin(angle)) * radius;
            vec2 v = (floor(uv * imageSize + offset) + 0.5) / imageSize;
            vec4 t = texture2D(surface, v);
            if (any(lessThan(v, vec2(0.0))) || any(greaterThan(v, vec2(1.0)))) t = vec4(0.0);
            bool boundary = (s.a < 0.0) != (t.a < 0.0);
            if (!boundary && outerOnly == 0 && s.a < 0.0 && t.a < 0.0) {
                vec3 delta = eyePosition(v, t.a) - p;
                // Tangent-plane continuity avoids outlining mesh triangulation.
                float step = abs(s.a - t.a);
                boundary = step > pixel * 12.0
                    && min(abs(dot(delta, s.xyz)), abs(dot(delta, t.xyz))) > pixel * 1.25;
                if (boundary && dot(s.xyz, t.xyz) >= 0.7 && step <= pixel * 250.0) {
                    vec2 direction = v - uv;
                    float before = texture2D(surface, uv - direction).a;
                    float after = texture2D(surface, v + direction).a;
                    // Parallel occluding surfaces jump beyond either local slope.
                    boundary = before < 0.0 && after < 0.0
                        && step > 4.0 * max(abs(s.a - before), abs(t.a - after)) + pixel;
                }
            }
            if (boundary) {
                ink = 1.0;
                if (t.a < 0.0 && (inkZ == 0.0 || t.a > inkZ)) {
                    inkZ = t.a;
                    z = min(z, texture2D(depth, v).r);
                }
            }
        }
    }
    if (ink > 0.0) {
        float fog = clamp((fogRange.y + inkZ) / max(fogRange.y - fogRange.x, 0.001), 0.0, 1.0);
        c = vec4(mix(background, edgeColor, fog), 1.0);
    }
    resultDepth = z;
    return c;
}
void main() {
    vec2 center = (gl_FragCoord.xy - viewportOrigin) * sampleScale;
    vec4 color = vec4(0.0);
    float nearest = 1.0;
    for (int y = -1; y <= 1; ++y) {
        for (int x = -1; x <= 1; ++x) {
            float z;
            color += sampleColor((center + vec2(float(x), float(y))) / imageSize, z) / 9.0;
            nearest = min(nearest, z);
        }
    }
    if (color.a == 0.0) discard;
    gl_FragColor = color;
    gl_FragDepth = nearest;
}
