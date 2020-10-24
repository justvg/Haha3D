#version 330 core
layout (location = 0) in vec3 aP;

out vec2 TexCoords;

void main()
{
    TexCoords = 0.5*(aP.xy + vec2(1.0, 1.0));
    gl_Position = vec4(aP, 1.0);
}