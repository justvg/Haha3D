#type VS
#version 330 core
layout (location = 0) in vec3 aP;
layout (location = 1) in vec3 aColor;

uniform mat4 Projection = mat4(1.0);
uniform mat4 View = mat4(1.0);
uniform mat4 Model = mat4(1.0);

out vec3 Color;

void main()
{
    Color = aColor;

    gl_Position = Projection * View * Model * vec4(aP, 1.0);
}

#type FS
#version 330 core
out vec4 FragColor;

in vec3 Color;

void main()
{
    FragColor = vec4(Color, 1.0);
}
