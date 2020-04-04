#include "haha3d_renderer_opengl.h"

shader::shader(char *VSSourceCode, char *FSSourceCode)
{
    GLuint VS = glCreateShader(GL_VERTEX_SHADER);
    GLuint FS = glCreateShader(GL_FRAGMENT_SHADER);
    ID = glCreateProgram();

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
}

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

shader::shader(char *ShaderFilepath)
{
    std::string SourceCode = ReadEntireFileAsString(ShaderFilepath);
    std::string VSSourceCode, FSSourceCode;

    char *TypeSpecifier = "#type";
    size_t VSTypeSpecifierP = SourceCode.find(TypeSpecifier, 0);
    size_t EndOfLine = SourceCode.find_first_of("\r\n", VSTypeSpecifierP);

    size_t VSSourceCodeStartP = SourceCode.find_first_not_of("\r\n", EndOfLine);
    size_t FSTypeSpecifierP = SourceCode.find(TypeSpecifier, VSSourceCodeStartP);
    EndOfLine = SourceCode.find_first_of("\r\n", FSTypeSpecifierP);
    size_t FSSourceCodeStartP = SourceCode.find_first_not_of("\r\n", EndOfLine);

    VSSourceCode = SourceCode.substr(VSSourceCodeStartP, FSTypeSpecifierP - VSSourceCodeStartP);
    FSSourceCode = SourceCode.substr(FSSourceCodeStartP, SourceCode.size() - FSSourceCodeStartP);

    *this = shader((char *)VSSourceCode.c_str(), (char *)FSSourceCode.c_str());
}

void shader::Use()
{
    glUseProgram(ID);
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
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(r32), (void *)0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
} 

PLATFORM_COMPILE_SHADER(CompileShader)
{
    shader Shader(ShaderFilepath);
    *Handle = (void *)((u64)Shader.ID);
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
                RenderCommandBuffer->CurrentShaderID = (GLuint)Command->Handle;
                glUseProgram(RenderCommandBuffer->CurrentShaderID);
            } break;
            
            case RenderCommand_PushMat4:
            {
                glUniformMatrix4fv(glGetUniformLocation(RenderCommandBuffer->CurrentShaderID, (char *)Command->Name), 1, GL_FALSE, Command->Elements);
            } break;
        }
    }
}