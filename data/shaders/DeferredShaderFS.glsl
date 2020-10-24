#version 330 core
out vec4 FragColor;

uniform sampler2D gPos;
uniform sampler2D gNormals;
uniform sampler2D gAlbedo;
uniform sampler2D AOTexture;

in vec2 TexCoords;

uniform vec3 LightDirView;

void main()
{
    vec3 Pos = texture(gPos, TexCoords).rgb;
    vec3 Normal = texture(gNormals, TexCoords).rgb;
    vec3 Albedo = texture(gAlbedo, TexCoords).rgb;
    float AO = texture(AOTexture, TexCoords).r;

    vec3 Ambient = 0.8*Albedo*AO;
    vec3 Diffuse = max(dot(Normal, -LightDirView), 0.0)*Albedo;
    Diffuse = vec3(0.0f);

    FragColor = vec4(Ambient + Diffuse, 1.0);
}