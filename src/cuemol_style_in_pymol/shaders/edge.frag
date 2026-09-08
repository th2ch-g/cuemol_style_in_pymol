#version 120
uniform vec3 edgeColor;
varying float visible;
void main() {
    if (visible < 0.5) discard;
    gl_FragColor = vec4(edgeColor, 1.0);
}
