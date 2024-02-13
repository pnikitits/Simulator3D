#version 120

uniform sampler2D p3d_Texture0;

// Input from vertex shader
varying vec2 texcoord;

// Output to the screen
void main() {
  vec4 color = texture2D(p3d_Texture0, texcoord);
  gl_FragColor = color.bgra;
}
