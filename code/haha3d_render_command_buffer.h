#pragma once

enum render_command_type
{
    RenderCommand_Clear,
    RenderCommand_DrawModel,
    RenderCommand_PushShader,
    RenderCommand_PopShader,
    RenderCommand_PushMat4,

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

#define MAX_RENDER_COMMAND_BUFFER_ENTRIES 1024
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
