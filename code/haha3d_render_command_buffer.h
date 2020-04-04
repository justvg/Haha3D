#pragma once

enum render_command_type
{
    RenderCommand_Clear,
    RenderCommand_DrawModel,
    RenderCommand_PushShader,
    RenderCommand_PushMat4,
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
            r32 *Elements;
        };
    };

    render_command_buffer_entry() {}
};

#define MAX_RENDER_COMMAND_BUFFER_ENTRIES 1024
struct render_command_buffer
{
    u32 CurrentShaderID;

    u32 CommandCount;
    render_command_buffer_entry Entries[MAX_RENDER_COMMAND_BUFFER_ENTRIES];

    render_command_buffer() {}
};
