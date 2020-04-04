#pragma once

class shader
{
    private:
        GLuint ID;

        friend PLATFORM_COMPILE_SHADER(CompileShader);

    public:
        shader() {}
        shader(char *VSSourceCode, char *FSSourceCode); 
        shader(char *ShaderFilepath); 

        void Use();
        void SetMat4(char *Name, mat4 &Mat);
};

void shader::SetMat4(char *Name, mat4 &Mat)
{
    glUniformMatrix4fv(glGetUniformLocation(ID, Name), 1, GL_FALSE, Mat.E);
}
