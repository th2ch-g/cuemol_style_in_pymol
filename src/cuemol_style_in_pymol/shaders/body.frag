#version 120
uniform int material;
uniform int principled;
uniform vec4 pbr;
uniform vec4 materialLighting;
uniform vec4 materialFinish;
uniform int perspective;
uniform vec3 background;
uniform vec2 fogRange;
uniform int pencilPreview;
uniform vec4 liveViewport;
uniform float liveScale;
varying vec3 positionEye;
varying vec3 normalEye;
varying vec3 baseColor;
float lambda(float a2, float cosine) {
    float c2 = max(cosine * cosine, 1e-12);
    return 0.5 * (sqrt(1.0 + a2 * (1.0 - c2) / c2) - 1.0);
}
#include pencil
vec3 pencilColor(vec3 n, vec3 view, float fog) {
    float flash = n.z > 0.0 ? clamp((n.z + 0.5) / 1.5, 0.0, 1.0) : 0.0;
    float ndl = dot(n, normalize(vec3(1.0)));
    float key = ndl > 0.0 ? clamp((ndl + 0.5) / 1.5, 0.0, 1.0) : 0.0;
    float tone = 0.05 + 0.85 * 1.302 * (0.6 * flash + 0.4 * key);
    tone *= 1.0 - pow(1.0 - max(dot(n, view), 0.0), 3.5) * (1.0 - 0.35 * clamp(tone, 0.0, 1.0));
    tone = pow(clamp((1.0 - (1.0 - tone) * fog) / 1.2, 0.0, 1.0), 2.4);
    tone = tone <= 0.0031308 ? 12.92 * tone : 1.055 * pow(tone, 1.0 / 2.4) - 0.055;
    tone = mix(tone, 1.0, smoothstep(0.81, 0.86, tone));
    vec3 paper = vec3(240.0, 236.0, 221.0) / 255.0;
    vec3 ink = baseColor * (0.4 + 0.6 * tone);
    vec3 luma = vec3(0.2126, 0.7152, 0.0722);
    float inkLuma = dot(ink, luma), target = dot(paper, luma) - 0.15;
    if (inkLuma > target) ink *= target / max(inkLuma, 0.0001);
    vec2 xy = 3.0 / max(liveScale, 1.0) * vec2(gl_FragCoord.x - liveViewport.x, liveViewport.w - (gl_FragCoord.y - liveViewport.y));
    paper *= 1.0 - pencilLayer(xy, vec2(0.5735764, 0.8191520), tone, 0.92, 0.0) * (1.0 - ink);
    paper *= 1.0 - pencilLayer(xy, vec2(0.8191520, -0.5735764), tone, 0.62, 1.0) * (1.0 - ink * 0.74);
    paper *= 1.0 - pencilLayer(xy, vec2(0.1736482, 0.9848078), tone, 0.34, 2.0) * (1.0 - ink * 0.38);
    return paper;
}
void main() {
    if (material < 0) { gl_FragData[0] = vec4(baseColor, 1.0); return; }
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
    if (pencilPreview != 0) result = pencilColor(n, view, fog);
    else result = max(mix(background, result, fog), 0.0);
    gl_FragData[0] = vec4(result, 1.0);
    gl_FragData[1] = vec4(n, positionEye.z);
}
