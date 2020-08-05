#include "haha3d_platform.h"
#include "haha3d_intrinsics.h"
#include "haha3d_shader.h"

struct model
{
    u32 VertexCount;
    u32 IndexCount;

    GLuint VAO, VBO, EBO;
};

internal void
InitModel(model *Model, u32 Size, r32 *Vertices, u32 VertexCount, u32 Stride)
{
    Model->VertexCount = VertexCount;

    glGenVertexArrays(1, &Model->VAO);
    glGenBuffers(1, &Model->VBO);
    glBindVertexArray(Model->VAO);
    glBindBuffer(GL_ARRAY_BUFFER, Model->VBO);
    glBufferData(GL_ARRAY_BUFFER, Size, Vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, Stride, (void *)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, Stride, (void *)(3*sizeof(r32)));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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
    shader ParticleShader;

    model Cube, Quad, Sphere;

    u32 GameObjectCount;
    game_object *Hero;
    game_object GameObjects[1024];

    u32 CollisionCount;
    collision_data Collisions[8000];

    camera Camera;

    b32 IsInitialized;
};