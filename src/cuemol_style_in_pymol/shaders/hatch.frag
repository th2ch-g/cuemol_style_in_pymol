#version 120
varying vec3 normalEye;
varying vec3 positionEye;
varying vec3 baseColor;
void main() {
    gl_FragData[0] = vec4(normalize(normalEye), positionEye.z);
    gl_FragData[1] = vec4(baseColor, 1.0);
}
