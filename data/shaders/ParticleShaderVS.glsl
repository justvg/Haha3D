#version 330 core
layout (location = 0) in vec3 aP;

uniform mat4 Projection = mat4(1.0);
uniform mat4 View = mat4(1.0);

uniform vec3 WorldP;
uniform vec3 CameraRight;
uniform vec3 CameraUp;
uniform vec2 Scale;

void main()
{
    vec3 P = WorldP + aP.x*CameraRight*Scale.x + aP.y*CameraUp*Scale.y;

    gl_Position = Projection * View * vec4(P, 1.0);
}