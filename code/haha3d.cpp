#include "haha3d_platform.h"
#include "haha3d_intrinsics.h"
#include "haha3d_math.cpp"

platform_api Platform;

struct shader
{
    void *Handle;
};

struct model
{
    u32 VertexCount;

    void *Handle;
};

internal void
InitModel(model *Model, u32 Size, r32 *Vertices, u32 VertexCount, u32 Stride)
{
    Model->VertexCount = VertexCount;

    Platform.InitBuffers(&Model->Handle, Size, Vertices, Stride);
}

#include "haha3d_render_command_buffer.cpp"

extern "C" 
GAME_UPDATE_AND_RENDER(GameUpdateAndRender)
{
    Platform = PlatformAPI;

    static shader Shader;
    static b32 TriDataInitialized = false;
    static model Triangle = {};
    static model Cube = {};
    if(!TriDataInitialized)
    {
        Platform.CompileShader(&Shader.Handle, "shaders/shader.glsl");

        r32 TriVertices[] = 
        {
            -0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f,
            0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.5f, 0.0f, 0.0f, 1.0f, 0.0f,
        };
        
        InitModel(&Triangle, sizeof(TriVertices), TriVertices, 3, 6*sizeof(r32));

        r32 CubeVertices[] = {
            // Back face
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            -0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            // Front face
            -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            -0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            // Left face
            -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            // Right face
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            // Bottom face
            -0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            // Top face
            -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f,  0.5f , 0.5f, 0.5f, 0.5f, 0.5f, 
            0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
        };

        InitModel(&Cube, sizeof(CubeVertices), CubeVertices, 36, 6*sizeof(r32));

        TriDataInitialized = true;
    }

    static r32 CameraPitch = 0.0f;
    static r32 CameraHead = 0.0f;
    r32 CameraRotationSensetivity = 0.1f;
    CameraPitch -= Input->MouseYDisplacement*CameraRotationSensetivity;
    CameraHead -= Input->MouseXDisplacement*CameraRotationSensetivity;

    CameraPitch = CameraPitch > 89.0f ? 89.0f : CameraPitch;
    CameraPitch = CameraPitch < -89.0f ? -89.0f : CameraPitch;

    r32 PitchRadians = Radians(CameraPitch);
    r32 HeadRadians = Radians(CameraHead);

    r32 CameraDistanceFromHero = 5.0f;
    r32 FloorDistanceFromHero = CameraDistanceFromHero * Cos(-PitchRadians);
    
    r32 XOffsetFromHero = FloorDistanceFromHero * Sin(HeadRadians);
    r32 YOffsetFromHero = CameraDistanceFromHero * Sin(-PitchRadians);
    r32 ZOffsetFromHero = FloorDistanceFromHero * Cos(HeadRadians);
    vec3 CameraOffsetFromHero = vec3(XOffsetFromHero, YOffsetFromHero, ZOffsetFromHero);

    static r32 HeroRotation = 0.0f;

    r32 Theta = Degrees(ATan2(-CameraOffsetFromHero.z, -CameraOffsetFromHero.x)) - 90.0f;
    if(Input->MoveForward.EndedDown)
    {
        HeroRotation = Theta;
    }
    if(Input->MoveBack.EndedDown)
    {
        HeroRotation = Theta + 180.0f;
    }
    if(Input->MoveRight.EndedDown)
    {
        HeroRotation = Theta - 90.0f;
    }
    if(Input->MoveLeft.EndedDown)
    {
        HeroRotation = Theta + 90.0f;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        HeroRotation = Theta - 45.0f;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        HeroRotation = Theta + 45.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        HeroRotation = Theta - 135.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        HeroRotation = Theta + 135.0f;
    }

    mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
    mat4 View = ViewRotationMatrixFromDirection(-CameraOffsetFromHero) * Translation(-CameraOffsetFromHero);
    mat4 Model = Rotation(HeroRotation, vec3(0.0f, 1.0f, 0.0f));

    Clear(RenderCommandBuffer, vec3(1.0f, 1.0f, 0.0f));

    PushShader(RenderCommandBuffer, Shader);
    PushMat4(RenderCommandBuffer, "Projection", &Projection);
    PushMat4(RenderCommandBuffer, "View", &View);
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &Cube);

    Model = Translation(vec3(0.0f, 0.0f, -3.0f));
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &Triangle);
}
