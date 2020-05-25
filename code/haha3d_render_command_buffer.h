#pragma once

enum render_command_type
{
    RenderCommand_Clear,
    RenderCommand_DrawModel,
    RenderCommand_DrawModelEBO,
    RenderCommand_PushShader,
    RenderCommand_PopShader,
    RenderCommand_PushMat4,
    RenderCommand_PushVec3,

    RenderCommand_EnableDepthTest,
    RenderCommand_DisableDepthTest,

    RenderCommand_DrawLine,
};

struct render_command_buffer_entry
{
    render_command_type Type;

    union
    {
        vec3 Color;
        struct
        {
            u32 VertexCount;
            void *Handle;
        };
        struct
        {
            u32 IndicesCount;
            void *Handle;
        };
        struct
        {
            char *Name;
            vec3 Vec3;
        };
        struct
        {
            char *Name;
            mat4 Matrix;
        };
        struct
        {
            vec3 From;
            vec3 To;
        };
    };

    render_command_buffer_entry() {}
};

#define MAX_RENDER_COMMAND_BUFFER_ENTRIES 4096
#define MAX_RENDER_COMMAND_BUFFER_SHADERS_IN_STACK 128
struct render_command_buffer
{
    u32 CurrentShaderID;
    u32 ShadersStackTop;
    u32 ShadersStack[MAX_RENDER_COMMAND_BUFFER_SHADERS_IN_STACK];

    u32 CommandCount;
    render_command_buffer_entry Entries[MAX_RENDER_COMMAND_BUFFER_ENTRIES];

    render_command_buffer() {}
};
