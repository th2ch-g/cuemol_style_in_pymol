#version 120
uniform int material;
varying vec3 normalEye;
varying vec3 positionWorld;
varying vec3 baseColor;
void main() {
    vec3 n = normalize(normalEye);
    vec3 light = normalize(vec3(0.35, 0.65, 1.0));
    float d = max(dot(n, light), 0.0);
    float s = max(dot(n, normalize(light + vec3(0.0, 0.0, 1.0))), 0.0);
    vec3 color = baseColor;
    float shade = 0.2 + 0.8 * d;
    if (material == 2) shade = 1.0;
    if (material == 1) shade = 0.75;
    if (material == 3) shade = 0.3 + 0.6 * d;
    if (material == 6 || material == 7) shade = 0.2 + 0.5 * d;
    if (material == 4) shade = d < 0.2 ? 0.48 : (d < 0.65 ? 0.76 : 1.0);
    if (material == 5) shade = d < 0.5 ? 0.5 : 1.0;
    if (material == 8 || material == 9) {
        shade = 0.3 + 0.7 * pow(0.5 + 0.5 * sin(n.y * 11.0 + n.z * 4.0), 2.0);
        color = (material == 9 ? vec3(1.0, 0.55, 0.27) : vec3(0.83, 0.9, 0.96)) * (0.65 + 0.35 * color);
    }
    if (material == 10) {
        float grain = sin(dot(positionWorld, vec3(12.1, 17.3, 8.7))) * sin(dot(positionWorld, vec3(29.7, 4.2, 21.3)));
        color *= 0.65 + 0.35 * grain;
    }
    if (material == 11 || material == 12) {
        float grain = 0.5 + 0.5 * sin(length(positionWorld.xz) * (material == 11 ? 3.0 : 7.0) + 0.5 * sin(positionWorld.y));
        color = vec3(0.25, 0.095, 0.035) + grain * vec3(0.55, 0.35, 0.15);
    }
    vec3 result = color * shade;
    if (material >= 6 && material <= 9) {
        result += 0.7 * pow(s, (material == 7 || material == 8) ? 70.0 : 25.0);
    }
    gl_FragColor = vec4(clamp(result, 0.0, 1.0), 1.0);
}
