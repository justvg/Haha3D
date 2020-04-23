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
    u32 IndexCount;

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
    vec3 ForceAccumulated;
    r32 TorqueAccumulated;

    // NOTE(georgy): Linear 
    r32 Mass;
    vec3 dP;
    vec3 P; // NOTE(georgy): This is center of mass in world coordinates

    // NOTE(georgy): Angular
    r32 MomentOfInertia;
    r32 Orientation; 
    r32 AngularSpeed;
};

enum game_object_type
{
    GameObject_Rectangle,
    GameObject_Cube,
};
struct game_object
{
    game_object_type Type;

    rigid_body RigidBody;
    model *Model;
    r32 Width, Height, Depth;
};

struct game_state
{
    shader Shader;

    model Cube, Quad, Sphere;

    u32 GameObjectCount;
    game_object *Hero;
    game_object GameObjects[8];

    camera Camera;

    b32 IsInitialized;
};
