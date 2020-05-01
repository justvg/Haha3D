#type VS
#version 330 core
layout (location = 0) in vec3 aP;
layout (location = 1) in vec3 aN;

uniform mat4 Projection = mat4(1.0);
uniform mat4 View = mat4(1.0);
uniform mat4 Model = mat4(1.0);

uniform vec3 CamP;

out vec3 FragPosWorld;
out vec3 CameraPosWorld;
out vec3 N;

void main()
{
    FragPosWorld = mat3(Model) * aP;
    CameraPosWorld = CamP;
    N = mat3(transpose(inverse(Model))) * aN;

    gl_Position = Projection * View * Model * vec4(aP, 1.0);
}

#type FS
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
