#version 330 core
out float AO;

uniform sampler2D AOTexture;

in vec2 TexCoords;

void main()
{
    vec2 TexelSize = 1.0 / vec2(textureSize(AOTexture, 0));

    float Result = 0.0f;
    for(int Y = -2; Y < 2; Y++)
    {
        for(int X = -2; X < 2; X++)
        {
            vec2 Offset = TexelSize*vec2(float(X), float(Y));
            Result += texture(AOTexture, TexCoords + Offset).r;
        }
    }

    AO = Result / (4.0 * 4.0);
}