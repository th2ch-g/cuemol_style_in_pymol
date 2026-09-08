#version 120
varying vec3 normalEye;
varying vec3 positionWorld;
varying vec3 baseColor;
void main() {
    vec4 eye = gl_ModelViewMatrix * gl_Vertex;
    gl_Position = gl_ProjectionMatrix * eye;
    gl_ClipVertex = eye;
    normalEye = normalize(gl_NormalMatrix * gl_Normal);
    positionWorld = gl_Vertex.xyz;
    baseColor = gl_Color.rgb;
}
