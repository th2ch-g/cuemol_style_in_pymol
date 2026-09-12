#version 120
uniform int material;
uniform int principled;
uniform vec4 pbr;
uniform vec4 materialLighting;
uniform vec4 materialFinish;
uniform int perspective;
uniform vec3 background;
uniform vec2 fogRange;
varying vec3 positionEye;
varying vec3 normalEye;
varying vec3 baseColor;
float lambda(float a2, float cosine) {
    float c2 = max(cosine * cosine, 1e-12);
    return 0.5 * (sqrt(1.0 + a2 * (1.0 - c2) / c2) - 1.0);
}
void main() {
    if (material < 0) { gl_FragColor = vec4(baseColor, 1.0); return; }
    vec3 n = normalize(normalEye);
    vec3 view = perspective == 0 ? vec3(0.0, 0.0, 1.0) : normalize(-positionEye);
    if (dot(n, view) < 0.0) n = -n;
    vec3 light = normalize(vec3(1.0));
    vec3 halfVector = normalize(light + view);
    float d = max(dot(n, light), 0.0);
    float flash = max(n.z, 0.0);
    float s = max(dot(n, halfVector), 0.0);
    vec3 result;
    if (principled != 0) {
        float alpha = max(pbr.y * pbr.y, 0.001);
        float a2 = alpha * alpha;
        float ndv = max(dot(n, view), 0.0001);
        vec3 f0 = mix(vec3(0.08 * pbr.z), baseColor, pbr.x);
        vec3 fresnel = f0 + (1.0 - f0) * pow(clamp(1.0 - dot(view, halfVector), 0.0, 1.0), 5.0);
        float t = s * s * (a2 - 1.0) + 1.0;
        float masking = 1.0 / (1.0 + lambda(a2, ndv) + lambda(a2, d));
        result = baseColor * (materialLighting.x + materialLighting.y * (1.0 - pbr.x) * (0.52 * d + 0.78 * flash));
        if (d > 0.0 && max(f0.r, max(f0.g, f0.b)) > 0.0)
            result += 0.52 * a2 / (t * t) * masking / (4.0 * ndv) * fresnel;
        result += (pbr.w > 0.0 ? vec3(pbr.w) : f0) * background;
    } else {
        float keyShape = d > 0.0 ? pow(d, materialFinish.x) : 0.0;
        float fillShape = flash > 0.0 ? pow(flash, materialFinish.x) : 0.0;
        result = baseColor * (materialLighting.x + materialLighting.y * (0.52 * keyShape + 0.78 * fillShape));
        if (d > 0.0) {
            float rv = max(dot(2.0 * d * n - light, view), 0.0);
            result += vec3(0.52 * materialFinish.y * pow(rv, materialFinish.z));
        }
    }
    float fog = clamp((fogRange.y + positionEye.z) / max(fogRange.y - fogRange.x, 0.001), 0.0, 1.0);
    gl_FragColor = vec4(max(mix(background, result, fog), 0.0), 1.0);
}
