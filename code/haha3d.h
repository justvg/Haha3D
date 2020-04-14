#include "haha3d_platform.h"
#include "haha3d_intrinsics.h"

global_variable platform_api Platform;

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

class camera
{
    public:
        r32 CameraPitch;
        r32 CameraHead;

        vec3 P;
        vec3 Dir;

        mat4 GetRotationMatrix(void);
};

struct rigid_body
{
    vec2 ForceAccumulated;
    r32 TorqueAccumulated;

    // NOTE(georgy): Linear 
    r32 Mass;
    vec2 dP;
    vec2 P; // NOTE(georgy): This is center of mass in world coordinates

    // NOTE(georgy): Angular
    r32 MomentOfInertia;
    r32 Orientation; 
    r32 AngularSpeed;
};
struct game_object
{
    rigid_body RigidBody;
};

struct game_state
{
    shader Shader;

    r32 HeroRotation;
    vec3 HeroP;
    model Cube, Quad;

    mat4 CubeOrientation;

    game_object HeroQuad;

    camera Camera;

    b32 IsInitialized;
};
