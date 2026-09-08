#version 120
varying vec3 normalEye;
varying vec3 positionEye;
varying vec3 baseColor;
void main() {
    vec4 eye = gl_ModelViewMatrix * gl_Vertex;
    positionEye = eye.xyz;
    normalEye = gl_NormalMatrix * gl_Normal;
    baseColor = gl_Color.rgb;
    gl_Position = gl_ProjectionMatrix * eye;
}
