#version 330 core
out float AO;

uniform sampler2D gPos;
uniform sampler2D gNormals;
uniform sampler2D gTangentNoise;

in vec2 TexCoords;

uniform mat4 Projection;

const float Radius = 1.0;
const int KernelSize = 32;
uniform vec3 SampleKernel[KernelSize];

uniform float WindowWidth;
uniform float WindowHeight;

void main()
{
    vec3 Pos = texture(gPos, TexCoords).rgb;
    vec3 Normal = texture(gNormals, TexCoords).rgb;

    vec2 NoiseTextureCoordsScale = vec2(WindowWidth, WindowHeight) / textureSize(gTangentNoise, 0);
    vec3 RandomVec = texture(gTangentNoise, NoiseTextureCoordsScale*TexCoords).rgb;

    vec3 Tangent = normalize(RandomVec - Normal*dot(RandomVec, Normal));
    vec3 Bitangent = cross(Normal, Tangent);
    mat3 TBN = mat3(Tangent, Bitangent, Normal);

    float Occlusion = 0.0f;
    for(int I = 0; I < KernelSize; I++)
    {
        vec3 SampleVec = TBN * SampleKernel[I];
        vec3 SampleP = Pos + SampleVec;

        vec4 Offset = Projection * vec4(SampleP, 1.0);
        Offset.xyz /= Offset.w;
        Offset.xyz = 0.5*Offset.xyz + vec3(0.5);
        float SampledDepth = texture(gPos, Offset.xy).z;

        float OcclusionFactor = smoothstep(0.0, 1.0, Radius / abs(Pos.z - SampledDepth));
        Occlusion += ((SampleP.z + 0.1 <= SampledDepth) ? 1.0 : 0.0) * OcclusionFactor;
    }

    AO = 1.0 - (Occlusion / KernelSize);
}