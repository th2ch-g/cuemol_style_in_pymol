#version 120
uniform sampler2D image;
uniform sampler2D depth;
uniform vec2 imageSize;
uniform vec2 viewportOrigin;
uniform int samples;
void main() {
    vec2 center = (gl_FragCoord.xy - viewportOrigin) * float(samples);
    vec4 color = vec4(0.0);
    float nearest = 1.0;
    if (samples == 1) {
        color = texture2D(image, center / imageSize);
        nearest = texture2D(depth, center / imageSize).r;
    } else {
    for (int y = -1; y <= 1; ++y) {
        for (int x = -1; x <= 1; ++x) {
            vec2 uv = (center + vec2(float(x), float(y))) / imageSize;
            color += texture2D(image, uv) / 9.0;
            nearest = min(nearest, texture2D(depth, uv).r);
        }
    }
    }
    if (color.a == 0.0) discard;
    gl_FragColor = color;
    gl_FragDepth = nearest;
}
