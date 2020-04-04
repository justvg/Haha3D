#include "haha3d_render_command_buffer.h"

inline b32
HasRoomForCommand(render_command_buffer *RenderCommandBuffer)
{
    b32 Result = (RenderCommandBuffer->CommandCount < MAX_RENDER_COMMAND_BUFFER_ENTRIES);
    return(Result);
}

internal render_command_buffer_entry *
PushCommand(render_command_buffer *RenderCommandBuffer, render_command_type Type)
{
    render_command_buffer_entry *Entry = 0;

    if(HasRoomForCommand(RenderCommandBuffer))
    {
        Entry = RenderCommandBuffer->Entries + RenderCommandBuffer->CommandCount++;
        Entry->Type = Type;
    }
    else
    {
        Assert(0);
    }

    return(Entry);
}

inline void
DrawModel(render_command_buffer *RenderCommandBuffer, model *Model)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_DrawModel);
    if(Entry)
    {
        Entry->VertexCount = Model->VertexCount;
        Entry->Handle = Model->Handle;
    }
}

inline void
Clear(render_command_buffer *RenderCommandBuffer, vec3 Color)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_Clear);
    if(Entry)
    {
        Entry->Color = Color;
    }
}

inline void
PushShader(render_command_buffer *RenderCommandBuffer, shader Shader)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_PushShader);
    if(Entry)
    {
        Entry->Handle = Shader.Handle;
    }
}

inline void
PushMat4(render_command_buffer *RenderCommandBuffer, char *Name, mat4 *Matrix)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_PushMat4);
    if(Entry)
    {
        Entry->Name = Name;
        Entry->Elements = Matrix->E;
    }
}