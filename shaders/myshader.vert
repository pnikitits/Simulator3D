#version 120

// Uniform inputs
uniform mat4 p3d_ModelViewProjectionMatrix;

// Vertex inputs
attribute vec4 p3d_Vertex;
attribute vec2 p3d_MultiTexCoord0;

// Output to fragment shader
varying vec2 texcoord;

void main() {
  gl_Position = p3d_ModelViewProjectionMatrix * p3d_Vertex;
  texcoord = p3d_MultiTexCoord0;
}
