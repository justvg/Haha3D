internal std::string
ReadEntireFileAsString(char *Filepath)
{
    std::string Result;

    std::ifstream File(Filepath, std::ios::in, std::ios::binary);
    if(File)
    {   
        File.seekg(0, std::ios::end);
        std::streamoff FileSize = File.tellg();
        File.seekg(0, std::ios::beg);

        Result.resize(FileSize);
        File.read(&Result[0], FileSize);
        File.close();
    }
    else
    {
        // TODO(georgy): Change OutputDebugString because it's Windows specific!
        OutputDebugString("Couldn't open file");
    }
    
    return(Result);
}

internal void 
InitOpenGLProperties(void)
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
}

internal
PLATFORM_INIT_BUFFERS(InitBuffers)
{
    GLuint *VAO = (GLuint *)Handle;
    GLuint *VBO = (GLuint *)((u8 *)Handle + sizeof(GLuint));

    glGenVertexArrays(1, VAO);
    glGenBuffers(1, VBO);
    glBindVertexArray(*VAO);
    glBindBuffer(GL_ARRAY_BUFFER, *VBO);
    glBufferData(GL_ARRAY_BUFFER, Size, Vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, Stride, (void *)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, Stride, (void *)(3*sizeof(r32)));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
} 

PLATFORM_COMPILE_SHADER(CompileShader)
{
    std::string SourceCode = ReadEntireFileAsString(ShaderFilepath);

    char *TypeSpecifier = "#type";
    size_t VSTypeSpecifierP = SourceCode.find(TypeSpecifier, 0);
    size_t EndOfLine = SourceCode.find_first_of("\r\n", VSTypeSpecifierP);

    size_t VSSourceCodeStartP = SourceCode.find_first_not_of("\r\n", EndOfLine);
    size_t FSTypeSpecifierP = SourceCode.find(TypeSpecifier, VSSourceCodeStartP);
    EndOfLine = SourceCode.find_first_of("\r\n", FSTypeSpecifierP);
    size_t FSSourceCodeStartP = SourceCode.find_first_not_of("\r\n", EndOfLine);

    std::string Str1, Str2;
    Str1 = SourceCode.substr(VSSourceCodeStartP, FSTypeSpecifierP - VSSourceCodeStartP);
    Str2 = SourceCode.substr(FSSourceCodeStartP, SourceCode.size() - FSSourceCodeStartP);
    char *VSSourceCode = (char *)Str1.c_str();
    char *FSSourceCode = (char *)Str2.c_str();

    GLuint VS = glCreateShader(GL_VERTEX_SHADER);
    GLuint FS = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint ID = glCreateProgram();

    GLint Success;
    char InfoLog[1024];

    glShaderSource(VS, 1, &VSSourceCode, 0);
    glCompileShader(VS);
    glGetShaderiv(VS, GL_COMPILE_STATUS, &Success);
    if(!Success)
    {
        glGetShaderInfoLog(VS, sizeof(InfoLog), 0, InfoLog);
        // TODO(georgy): Change OutputDebugString because it's Windows specific! Renderer is the third tier!
        OutputDebugString("ERROR::SHADER_COMPILATION_ERROR of type: VS\n");
        OutputDebugString(InfoLog);
        OutputDebugString("\n");
    }

    glShaderSource(FS, 1, &FSSourceCode, 0);
    glCompileShader(FS);
    glGetShaderiv(FS, GL_COMPILE_STATUS, &Success);
    if(!Success)
    {
        glGetShaderInfoLog(FS, sizeof(InfoLog), 0, InfoLog);
        // TODO(georgy): Change OutputDebugString because it's Windows specific! Renderer is the third tier!
        OutputDebugString("ERROR::SHADER_COMPILATION_ERROR of type: FS\n");
        OutputDebugString(InfoLog);
        OutputDebugString("\n");
    }

    glAttachShader(ID, VS);
    glAttachShader(ID, FS);
    glLinkProgram(ID);
    glGetProgramiv(ID, GL_LINK_STATUS, &Success);
    if(!Success)
    {
        glGetProgramInfoLog(ID, sizeof(InfoLog), 0, InfoLog);
        // TODO(georgy): Change OutputDebugString because it's Windows specific! Renderer is the third tier!
        OutputDebugString("ERROR::PROGRAM_LINKING_ERROR of type:: PROGRAM\n");
		OutputDebugString(InfoLog);
		OutputDebugString("\n");
    }

    glDeleteShader(VS);
    glDeleteShader(FS);

    *Handle = (void *)((u64)ID);
}

internal void
RenderCommands(render_command_buffer *RenderCommandBuffer)
{
    for(u32 CommandIndex = 0;
        CommandIndex < RenderCommandBuffer->CommandCount;
        CommandIndex++)
    {
        render_command_buffer_entry *Command = RenderCommandBuffer->Entries + CommandIndex;

        switch(Command->Type)
        {
            case RenderCommand_Clear:
            {
                glClearColor(Command->Color.x, Command->Color.y, Command->Color.z, 1.0f);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            } break;

            case RenderCommand_DrawModel:
            {
                GLuint VAO = (GLuint)Command->Handle;
                
                glBindVertexArray(VAO);
                glDrawArrays(GL_TRIANGLES, 0, Command->VertexCount);
                glBindVertexArray(0);
            } break;

            case RenderCommand_PushShader:
            {
                Assert(RenderCommandBuffer->ShadersStackTop < MAX_RENDER_COMMAND_BUFFER_SHADERS_IN_STACK);

                RenderCommandBuffer->CurrentShaderID = (GLuint)Command->Handle;
                RenderCommandBuffer->ShadersStack[RenderCommandBuffer->ShadersStackTop++] = RenderCommandBuffer->CurrentShaderID;
                glUseProgram(RenderCommandBuffer->CurrentShaderID);
            } break;

            case RenderCommand_PopShader:
            {
                Assert(RenderCommandBuffer->ShadersStackTop > 0);

                RenderCommandBuffer->ShadersStackTop--;
                 
                // TODO(georgy): I can get rid of this cmp, if I always push 0 in the 0-th slot of shader stack. Do I want to do this?
                if(RenderCommandBuffer->ShadersStackTop == 0)
                {
                    RenderCommandBuffer->CurrentShaderID = 0;    
                }
                else
                {
                    RenderCommandBuffer->CurrentShaderID = RenderCommandBuffer->ShadersStack[RenderCommandBuffer->ShadersStackTop - 1];
                }

                glUseProgram(RenderCommandBuffer->CurrentShaderID);
            } break;
            
            case RenderCommand_PushMat4:
            {
                glUniformMatrix4fv(glGetUniformLocation(RenderCommandBuffer->CurrentShaderID, (char *)Command->Name), 1, GL_FALSE, Command->Matrix.E);
            } break;
        }
    }
}