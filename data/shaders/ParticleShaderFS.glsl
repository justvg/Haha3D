#version 330 core
out vec4 FragColor;

uniform vec4 ColorMult;

void main()
{
    FragColor = ColorMult * vec4(1.0, 1.0, 1.0, 1.0);
}
