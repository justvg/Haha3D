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
#if 0
inline void
DrawModelEBO(render_command_buffer *RenderCommandBuffer, model *Model)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_DrawModelEBO);
    if(Entry)
    {
        Entry->IndicesCount = Model->IndexCount;
        Entry->Handle = Model->Handle;
    }
}
#endif
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
PopShader(render_command_buffer *RenderCommandBuffer)
{
    PushCommand(RenderCommandBuffer, RenderCommand_PopShader);
}

inline void
PushVec3(render_command_buffer *RenderCommandBuffer, char *Name, vec3 Vector)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_PushVec3);
    if(Entry)
    {
        Entry->Name = Name;
        Entry->Vec3 = Vector;
    }
}

inline void
PushMat4(render_command_buffer *RenderCommandBuffer, char *Name, mat4 *Matrix)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_PushMat4);
    if(Entry)
    {
        Entry->Name = Name;
        Entry->Matrix = *Matrix;
    }
}

inline void
DrawLine(render_command_buffer *RenderCommandBuffer, vec3 From, vec3 To)
{
    render_command_buffer_entry *Entry = PushCommand(RenderCommandBuffer, RenderCommand_DrawLine);
    if(Entry)
    {
        Entry->From = From;
        Entry->To = To;
    }
}