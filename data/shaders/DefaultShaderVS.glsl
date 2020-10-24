#version 330 core
layout (location = 0) in vec3 aP;
layout (location = 1) in vec3 aN;

uniform mat4 Projection = mat4(1.0);
uniform mat4 View = mat4(1.0);
uniform mat4 Model = mat4(1.0);

out VS_OUT
{
    vec3 FragPosView;
    vec3 NormalView;
} vs_out;

void main()
{
    mat4 ViewModel = View * Model;
    vs_out.FragPosView = vec3(ViewModel * vec4(aP, 1.0));
    vs_out.NormalView = mat3(transpose(inverse(ViewModel))) * aN;

    gl_Position = Projection * View * Model * vec4(aP, 1.0);
}