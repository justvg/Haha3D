
internal void
GameUpdateAndRender(game_input *Input, u32 WindowWidth, u32 WindowHeight)
{
    static shader Shader("shaders/Shader.glsl");

    static b32 TriDataInitialized = false;
    static GLuint TriVAO = 0, TriVBO = 0;
    static GLuint CubeVAO = 0, CubeVBO = 0;
    if(!TriDataInitialized)
    {
        r32 TriVertices[] = 
        {
            -0.5f, -0.5f, 0.0f,
            0.5f, -0.5f, 0.0f,
            0.0f, 0.5f, 0.0f
        };

        glGenVertexArrays(1, &TriVAO);
        glGenBuffers(1, &TriVBO);
        glBindVertexArray(TriVAO);
        glBindBuffer(GL_ARRAY_BUFFER, TriVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(TriVertices), TriVertices, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(r32), (void *)0);
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        
        r32 CubeVertices[] = {
            // Back face
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            -0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            // Front face
            -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            -0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            // Left face
            -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            // Right face
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            // Bottom face
            -0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            // Top face
            -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f,  0.5f , 0.5f, 0.5f, 0.5f, 0.5f, 
            0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
        };
        glGenVertexArrays(1, &CubeVAO);
        glGenBuffers(1, &CubeVBO);
        glBindVertexArray(CubeVAO);
        glBindBuffer(GL_ARRAY_BUFFER, CubeVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(CubeVertices), CubeVertices, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(r32), (void *)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(r32), (void *)(3*sizeof(r32)));
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        TriDataInitialized = true;
    }

    static r32 CameraPitch = 0.0f;
    static r32 CameraHead = 0.0f;
    r32 CameraRotationSensetivity = 0.1f;
    CameraPitch -= Input->MouseYDisplacement*CameraRotationSensetivity;
    CameraHead -= Input->MouseXDisplacement*CameraRotationSensetivity;

    CameraPitch = CameraPitch > 89.0f ? 89.0f : CameraPitch;
    CameraPitch = CameraPitch < -89.0f ? -89.0f : CameraPitch;

    r32 PitchRadians = Radians(CameraPitch);
    r32 HeadRadians = Radians(CameraHead);

    r32 CameraDistanceFromHero = 5.0f;
    r32 FloorDistanceFromHero = CameraDistanceFromHero * Cos(-PitchRadians);
    
    r32 XOffsetFromHero = FloorDistanceFromHero * Sin(HeadRadians);
    r32 YOffsetFromHero = CameraDistanceFromHero * Sin(-PitchRadians);
    r32 ZOffsetFromHero = FloorDistanceFromHero * Cos(HeadRadians);
    vec3 CameraOffsetFromHero = vec3(XOffsetFromHero, YOffsetFromHero, ZOffsetFromHero);

    static r32 HeroRotation = 0.0f;

    r32 Theta = Degrees(ATan2(-CameraOffsetFromHero.z, -CameraOffsetFromHero.x)) - 90.0f;
    if(Input->MoveForward.EndedDown)
    {
        HeroRotation = Theta;
    }
    if(Input->MoveBack.EndedDown)
    {
        HeroRotation = Theta + 180.0f;
    }
    if(Input->MoveRight.EndedDown)
    {
        HeroRotation = Theta - 90.0f;
    }
    if(Input->MoveLeft.EndedDown)
    {
        HeroRotation = Theta + 90.0f;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        HeroRotation = Theta - 45.0f;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        HeroRotation = Theta + 45.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        HeroRotation = Theta - 135.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        HeroRotation = Theta + 135.0f;
    }

    mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
    mat4 View = ViewRotationMatrixFromDirection(-CameraOffsetFromHero) * Translation(-CameraOffsetFromHero);
    mat4 Model = Rotation(HeroRotation, vec3(0.0f, 1.0f, 0.0f));

    Shader.Use();
    Shader.SetMat4("Projection", Projection);
    Shader.SetMat4("View", View);
    Shader.SetMat4("Model", Model);
    glBindVertexArray(CubeVAO);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
    Model = Translation(vec3(0.0f, 0.0f, -3.0f));
    Shader.SetMat4("Model", Model);
    glBindVertexArray(TriVAO);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
}
