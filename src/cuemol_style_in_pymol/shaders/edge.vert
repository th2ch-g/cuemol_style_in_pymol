#version 120
uniform int creases;
uniform float edgeWidth;
uniform vec2 viewportSize;
varying float visible;
void main() {
    vec4 eye = gl_ModelViewMatrix * gl_Vertex;
    vec3 a = normalize(gl_NormalMatrix * gl_Normal);
    vec3 b = normalize(gl_NormalMatrix * gl_MultiTexCoord0.xyz);
    vec3 viewDir = gl_ProjectionMatrix[3][3] == 0.0 ? normalize(-eye.xyz) : vec3(0.0, 0.0, 1.0);
    float fa = dot(a, viewDir);
    float fb = dot(b, viewDir);
    visible = (fa * fb <= 0.0 || (creases == 1 && dot(a, b) < 0.5 && max(fa, fb) > 0.0)) ? 1.0 : 0.0;
    gl_Position = gl_ProjectionMatrix * eye;
    vec4 otherEye = gl_ModelViewMatrix * vec4(gl_MultiTexCoord1.xyz, 1.0);
    vec4 otherClip = gl_ProjectionMatrix * otherEye;
    vec2 delta = (otherClip.xy / otherClip.w - gl_Position.xy / gl_Position.w) * viewportSize;
    vec2 side = vec2(-delta.y, delta.x) / max(length(delta), 0.00001);
    float distance = gl_ProjectionMatrix[3][3] == 0.0 ? max(abs(eye.z), 0.00001) : 1.0;
    float pixels = max(1.0, edgeWidth * abs(gl_ProjectionMatrix[1][1]) * viewportSize.y / (2.0 * distance));
    gl_Position.xy += side * gl_MultiTexCoord1.w * pixels / viewportSize * gl_Position.w;
    gl_ClipVertex = eye;
}
