#version 330 core
out vec4 FragColor;

uniform vec3 Color;

in vec3 FragPosWorld;
in vec3 CameraPosWorld;
in vec3 N;

void main()
{
    const vec3 LightDir = normalize(vec3(-0.5, -0.7, -0.5));
    vec3 Normal = normalize(N);

    vec3 Ambient = 0.2 * Color;
    vec3 Diffuse = max(dot(Normal, -LightDir), 0.0) * Color;

    FragColor = vec4(Ambient + Diffuse, 1.0);
}