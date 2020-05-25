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
    vec3 TorqueAccumulated;

    // NOTE(georgy): Linear 
    r32 Mass;
    vec3 dP;
    vec3 P; // NOTE(georgy): This is center of mass in world coordinates

    // NOTE(georgy): Angular
    mat3 InertiaTensor;
    mat3 InverseInertiaTensor;
    mat3 GlobalInverseInertiaTensor;
    mat3 Orientation; 
    mat3 InverseOrientation;
    vec3 AngularMomentum;
    vec3 AngularVelocity;

    r32 CoeffOfRestitution;
    r32 CoeffOfFriction;
};

enum game_object_type
{
    GameObject_Rectangle,
    GameObject_Cube,
    GameObject_Wall,
};
struct game_object
{
    game_object_type Type;

    rigid_body RigidBody;
    model *Model;
    r32 Width, Height, Depth;
};

struct contact_point
{
    b32 Persistent;

    vec3 GlobalA;
    vec3 LocalA;
    vec3 GlobalB;
    vec3 LocalB;
    vec3 Normal;
    r32 Penetration;

    r32 AccumulatedImpulse;
    r32 AccumulatedImpulseTangent1;
    r32 AccumulatedImpulseTangent2;
};
struct collision_manifold
{
    u32 ContactCount;
    // NOTE(georgy): We use only 4. But we need 1 extra slot at the end
    contact_point Contacts[5];
};

struct collision_data
{
    b32 CollisionThisFrame;
    game_object *A;
    game_object *B;
    collision_manifold Manifold;
};

struct game_state
{
    shader Shader;
    shader DebugShader;

    model Cube, Quad, Sphere;

    u32 GameObjectCount;
    game_object *Hero;
    game_object GameObjects[1024];

    u32 CollisionCount;
    collision_data Collisions[8000];

    camera Camera;

    b32 IsInitialized;
};