#include "haha3d.h"
#include "haha3d_math.cpp"
#include "haha3d_render_command_buffer.cpp"

mat4 
camera::GetRotationMatrix(void)
{ 
    mat4 Result = LookAt(P, P + Dir); 

    return(Result);
} 

internal vec2
GetPointOfContact(rigid_body *RigidBody, r32 t, vec2 DeltaP, 
                  r32 HalfWidth, vec2 OrientationX, r32 HalfHeight, vec2 OrientationY,
                  plane Plane, r32 Radius)
{
    vec2 RigidBodyCorners[4];
    RigidBodyCorners[0] = RigidBody->P + t*DeltaP + 
        (HalfWidth*OrientationX + HalfHeight*OrientationY);
    RigidBodyCorners[1] = RigidBody->P + t*DeltaP + 
        (HalfWidth*OrientationX - HalfHeight*OrientationY);
    RigidBodyCorners[2] = RigidBody->P + t*DeltaP + 
        (-HalfWidth*OrientationX - HalfHeight*OrientationY);
    RigidBodyCorners[3] = RigidBody->P + t*DeltaP + 
        (-HalfWidth*OrientationX + HalfHeight*OrientationY);

    u32 ClosestCount = 0;
    u32 ClosestIndecies[4];

    u32 TheClosestToPlaneIndex = 0;
    r32 ClosestDistance = 1000000.0f;

    for(u32 Corner = 0;
        Corner < 4;
        Corner++)
    {
        r32 Distance = Absolute(Dot(RigidBodyCorners[Corner], Plane.N) - Plane.D);
        if(Distance > 0.0f)
        {
            r32 Diff = ClosestDistance - Distance;
            if(Absolute(Diff) <= Radius*0.1f)
            {
                ClosestIndecies[ClosestCount++] = Corner;
            }
            else if(Diff > 0.0f)
            {
                TheClosestToPlaneIndex = Corner;
                ClosestDistance = Distance;

                ClosestCount = 0;
                ClosestIndecies[ClosestCount++] = Corner;
            }
        }
    }

    // Assert(ClosestCount <= 2);

    vec2 PointOfContact;
    if(ClosestCount > 2)
    {
        PointOfContact = Lerp(RigidBodyCorners[ClosestIndecies[0]], RigidBodyCorners[ClosestIndecies[1]], 0.5f);
    }
    else
    {
        PointOfContact = RigidBodyCorners[TheClosestToPlaneIndex];
    }

    return(PointOfContact);
}

internal void
ComputeProjectionInterval(vec2 *Vertices, u32 Count, vec2 V, r32 &Min, r32 &Max)
{
    Min = Max = Dot(Vertices[0], V);
    for(u32 Vertex = 1;
        Vertex < Count;
        Vertex++)
    {
        r32 Value = Dot(Vertices[Vertex], V);
        if(Value < Min)
        {
            Min = Value;
        }
        else if(Value > Max)
        {
            Max = Value;
        }
    }
}

internal b32 
NoIntersection(r32 MinA, r32 MaxA, r32 MinB, r32 MaxB, r32 Speed, r32 &tFirst, r32 &tLast)
{
    if(MaxA < MinB)
    {
        // NOTE(georgy): Interval of A is initially on 'left' of interval B
        
        if(Speed <= 0.0f) { return(true); }

        r32 t = (MinB - MaxA) / Speed;
        if(t > tFirst)
        {
            tFirst = t;
        }

        t = (MaxB - MinA) / Speed;
        if(t < tLast)
        {
            tLast = t;
        }
    }
    else if(MaxB < MinA)
    {
        // NOTE(georgy): Interval of A is initially on 'right' of interval B

        if(Speed >= 0) { return(true); }

        r32 t = (MaxB - MinA) / Speed;
        if(t > tFirst)
        {
            tFirst = t;
        }

        t = (MinB - MaxA) / Speed;
        if(t < tLast)
        {
            tLast = t;
        }
    }
    else
    {
        // NOTE(georgy): Interval of A and interval of B initially overlap

        if(Speed > 0)
        {
            r32 t = (MaxB - MinA) / Speed;
            if(t < tLast)
            {
                tLast = t;
            }
        }
        else if(Speed < 0)
        {
            r32 t = (MinB - MaxA) / Speed;
            if(t < tLast)
            {
                tLast = t;
            }
        }
    }

    b32 NoIntersect = (tFirst > tLast);
    return(NoIntersect);
}

internal b32
TestIntersection(game_object *A, vec2 DeltaP, game_object *B, r32 *t)
{
    // NOTE(georgy): Process as if B is stationary, A is moving

    r32 tFirst = 0.0f;
    r32 tLast = 10000.0f;

    vec2 AAxisX = vec2(Cos(Radians(A->RigidBody.Orientation)), Sin(Radians(A->RigidBody.Orientation)));
    vec2 AAxisY = Perp(AAxisX);
    // TODO(georgy): This function works for any convex polygon!
    // NOTE(georgy): CCW order
    vec2 APolygonVertices[4] = 
    {
        vec2(A->RigidBody.P + 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY),
        vec2(A->RigidBody.P - 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY),
        vec2(A->RigidBody.P - 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY),
        vec2(A->RigidBody.P + 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY),
    };

    vec2 BAxisX = vec2(Cos(Radians(B->RigidBody.Orientation)), Sin(Radians(B->RigidBody.Orientation)));
    vec2 BAxisY = Perp(BAxisX);
    // TODO(georgy): This function works for any convex polygon!
    // NOTE(georgy): CCW order
    vec2 BPolygonVertices[4] = 
    {
        vec2(B->RigidBody.P + 0.5f*B->Width*BAxisX + 0.5f*B->Height*BAxisY),
        vec2(B->RigidBody.P - 0.5f*B->Width*BAxisX + 0.5f*B->Height*BAxisY),
        vec2(B->RigidBody.P - 0.5f*B->Width*BAxisX - 0.5f*B->Height*BAxisY),
        vec2(B->RigidBody.P + 0.5f*B->Width*BAxisX - 0.5f*B->Height*BAxisY),
    };

    // NOTE(georgy): Test normals of A's edges for separation
    for(u32 I0 = 0, I1 = ArrayCount(APolygonVertices) - 1;
        I0 < ArrayCount(APolygonVertices);
        I1 = I0, I0++)
    {
        vec2 Edge = APolygonVertices[I0] - APolygonVertices[I1];
        vec2 Normal = vec2(Edge.y, -Edge.x);
        r32 Speed = Dot(DeltaP, Normal);

        r32 MinA, MaxA, MinB, MaxB;
        ComputeProjectionInterval(APolygonVertices, ArrayCount(APolygonVertices), Normal, MinA, MaxA);
        ComputeProjectionInterval(BPolygonVertices, ArrayCount(BPolygonVertices), Normal, MinB, MaxB);

        if(NoIntersection(MinA, MaxA, MinB, MaxB, Speed, tFirst, tLast))
        {
            return(false);
        }
    }

    // NOTE(georgy): Test normals of B's edges for separation
    for(u32 I0 = 0, I1 = ArrayCount(BPolygonVertices) - 1;
        I0 < ArrayCount(BPolygonVertices);
        I1 = I0, I0++)
    {
        vec2 Edge = BPolygonVertices[I0] - BPolygonVertices[I1];
        vec2 Normal = vec2(Edge.y, -Edge.x);
        r32 Speed = Dot(DeltaP, Normal);

        r32 MinA, MaxA, MinB, MaxB;
        ComputeProjectionInterval(APolygonVertices, ArrayCount(APolygonVertices), Normal, MinA, MaxA);
        ComputeProjectionInterval(BPolygonVertices, ArrayCount(BPolygonVertices), Normal, MinB, MaxB);

        if(NoIntersection(MinA, MaxA, MinB, MaxB, Speed, tFirst, tLast))
        {
            return(false);
        }
    }

    *t = tFirst;
    return(true);
}

extern "C" 
GAME_UPDATE_AND_RENDER(GameUpdateAndRender)
{
    Platform = Memory->PlatformAPI;

    Assert(sizeof(game_state) <= Memory->PermanentStorageSize);
    game_state *GameState = (game_state *)Memory->PermanentStorage;
    if(!GameState->IsInitialized)
    {
        Platform.CompileShader(&GameState->Shader.Handle, "shaders/shader.glsl");

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

        InitModel(&GameState->Cube, sizeof(CubeVertices), CubeVertices, 36, 6*sizeof(r32));

        r32 QuadVertices[] = 
        {
            -0.5f, 0.5f, 0.0f, 1.0f, 1.0f, 1.0f,
            -0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 1.0f,
            0.5f, 0.5f, 0.0f, 1.0f, 1.0f, 1.0f,

            0.5f, 0.5f, 0.0f, 1.0f, 1.0f, 1.0f,
            -0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 1.0f,
            0.5f, -0.5f, 0.0f, 1.0f, 1.0f, 1.0f,
        };

        InitModel(&GameState->Quad, sizeof(QuadVertices), QuadVertices, 6, 6*sizeof(r32));

        GameState->CubeOrientation = Identity();

        GameState->GameObjectCount = 2;
        GameState->Hero = GameState->GameObjects;
        GameState->Hero->Type = GameObject_Rectangle;
        GameState->Hero->Model = &GameState->Quad;
        GameState->Hero->Width = 0.1f;
        GameState->Hero->Height = 0.05f;
        GameState->Hero->RigidBody.Mass = 3.0f;
        GameState->Hero->RigidBody.MomentOfInertia = (1.0f / 12.0f) * GameState->Hero->RigidBody.Mass * 
                                                      (Square(GameState->Hero->Width) + Square(GameState->Hero->Height));

        // GameState->GameObjects[1].Type = GameObject_Rectangle;
        // GameState->GameObjects[1].Model = &GameState->Quad;
        // GameState->GameObjects[1].RigidBody.P = vec2(0.2f, 0.0f);
        // GameState->GameObjects[1].Width = 0.1f;
        // GameState->GameObjects[1].Height = 0.05f;
        // GameState->GameObjects[1].RigidBody.Mass = 3.0f;
        // GameState->GameObjects[1].RigidBody.MomentOfInertia = (1.0f / 12.0f) * GameState->Hero->RigidBody.Mass * 
        //                                                       (Square(GameState->GameObjects[1].Width) + Square(GameState->GameObjects[1].Height));

        GameState->GameObjects[1].Type = GameObject_Rectangle;
        GameState->GameObjects[1].Model = &GameState->Quad;
        GameState->GameObjects[1].RigidBody.P = vec2(0.0f, -0.35f);
        GameState->GameObjects[1].Width = 0.1f;
        GameState->GameObjects[1].Height = 0.1f;
        GameState->GameObjects[1].RigidBody.Mass = 5.0f;
        GameState->GameObjects[1].RigidBody.MomentOfInertia = (1.0f / 12.0f) * GameState->Hero->RigidBody.Mass * 
                                                              (Square(GameState->GameObjects[1].Width) + Square(GameState->GameObjects[1].Height));

        GameState->IsInitialized = true;
    }

    r32 dt = Input->dtForFrame;

    camera *Camera = &GameState->Camera;
#if 0
    r32 CameraRotationSensetivity = 0.1f;
    Camera->CameraPitch -= Input->MouseYDisplacement*CameraRotationSensetivity;
    Camera->CameraHead -= Input->MouseXDisplacement*CameraRotationSensetivity;

    Camera->CameraPitch = Camera->CameraPitch > 89.0f ? 89.0f : Camera->CameraPitch;
    Camera->CameraPitch = Camera->CameraPitch < -89.0f ? -89.0f : Camera->CameraPitch;

    r32 PitchRadians = Radians(Camera->CameraPitch);
    r32 HeadRadians = Radians(Camera->CameraHead);

    r32 CameraDistanceFromHero = 5.0f;
    r32 FloorDistanceFromHero = CameraDistanceFromHero * Cos(-PitchRadians);
    
    r32 XOffsetFromHero = FloorDistanceFromHero * Sin(HeadRadians);
    r32 YOffsetFromHero = CameraDistanceFromHero * Sin(-PitchRadians);
    r32 ZOffsetFromHero = FloorDistanceFromHero * Cos(HeadRadians);
    vec3 CameraOffsetFromHero = vec3(XOffsetFromHero, YOffsetFromHero, ZOffsetFromHero);

    vec3 CameraForward = Normalize(-CameraOffsetFromHero);
    vec3 CameraRight = Normalize(Cross(CameraForward, vec3(0.0f, 1.0f, 1.0f)));
    r32 Theta = Degrees(ATan2(CameraForward.z, CameraForward.x)) - 90.0f;
    if(Input->MoveForward.EndedDown)
    {
        GameState->HeroRotation = Theta;
        GameState->HeroP += 20.0f*dt*CameraForward;
    }
    if(Input->MoveBack.EndedDown)
    {
        GameState->HeroRotation = Theta + 180.0f;
        GameState->HeroP += 20.0f*dt*-CameraForward;
    }
    if(Input->MoveRight.EndedDown)
    {
        GameState->HeroRotation = Theta - 90.0f;
        GameState->HeroP += 20.0f*dt*CameraRight;
    }
    if(Input->MoveLeft.EndedDown)
    {
        GameState->HeroRotation = Theta + 90.0f;
        GameState->HeroP += 20.0f*dt*-CameraRight;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        GameState->HeroRotation = Theta - 45.0f;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        GameState->HeroRotation = Theta + 45.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        GameState->HeroRotation = Theta - 135.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        GameState->HeroRotation = Theta + 135.0f;
    }

    Camera->P = GameState->HeroP + CameraOffsetFromHero;
    Camera->Dir = CameraForward;
#else

    GameState->Hero->RigidBody.ForceAccumulated = vec2(0.0f, 0.0f);
    GameState->Hero->RigidBody.TorqueAccumulated = 0.0f;

    static b32 ForceActsOnCM = true;
    if(WasDown(&Input->MouseLeft))
    {
        ForceActsOnCM = !ForceActsOnCM;
    }

    static b32 GravityEnabled = false;
    if(WasDown(&Input->MouseRight))
    {
        GravityEnabled = !GravityEnabled;
    }

    rigid_body *HeroRigidBody = &GameState->Hero->RigidBody;
    r32 HeroHalfWidth = 0.5f*GameState->Hero->Width;
    r32 HeroHalfHeight = 0.5f*GameState->Hero->Height;
    vec2 HeroOrientationAxisX = vec2(Cos(Radians(HeroRigidBody->Orientation)), Sin(Radians(HeroRigidBody->Orientation)));
    vec2 HeroOrientationAxisY = Perp(HeroOrientationAxisX);

    // NOTE(georgy): At the moment these forces are acting only on the CM of the hero
    if(Input->MoveForward.EndedDown)
    {
        vec2 Force = 200.0f*vec2(0.0f, 1.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(HeroHalfHeight*HeroOrientationAxisY + 0.2f*HeroHalfWidth*HeroOrientationAxisX, Force);
        }
    }
    if(Input->MoveBack.EndedDown)
    {
        vec2 Force = -200.0f*vec2(0.0f, 1.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(-HeroHalfHeight*HeroOrientationAxisY - 0.2f*HeroHalfWidth*HeroOrientationAxisX, Force);
        }
    }
    if(Input->MoveRight.EndedDown)
    {
        vec2 Force = 200.0f*vec2(1.0f, 0.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(-HeroHalfWidth*HeroOrientationAxisX + 0.2f*HeroHalfHeight*HeroOrientationAxisY, Force);
        }
    }
    if(Input->MoveLeft.EndedDown)
    {
        vec2 Force = -200.0f*vec2(1.0f, 0.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(HeroHalfWidth*HeroOrientationAxisX - 0.2f*HeroHalfHeight*HeroOrientationAxisY, Force);
        }
    }

    Camera->P = vec3(0.0f, 0.0f, 1.0f);
    Camera->Dir = vec3(0.0f, 0.0f, -1.0f);

#endif
    mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
    mat4 View = Camera->GetRotationMatrix();

    Clear(RenderCommandBuffer, vec3(1.0f, 0.0f, 0.0f));

    PushShader(RenderCommandBuffer, GameState->Shader);
    PushMat4(RenderCommandBuffer, "Projection", &Projection);
    PushMat4(RenderCommandBuffer, "View", &View);
#if 0
    mat4 Model = Translation(GameState->HeroP) * Rotation(GameState->HeroRotation, vec3(0.0f, 1.0f, 0.0f));
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Cube);

    Model = Translation(vec3(0.0f, -2.0f, 0.0f)) * Rotation(-90.0f, vec3(1.0f, 0.0f, 0.0f)) * Scaling(6.0f);
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Plane);

    GameState->CubeOrientation = Rotation(dt*180.0f, vec3(0.0f, 1.0f, 0.0f)) * GameState->CubeOrientation;
    Model = Translation(vec3(0.0f, -1.5f, 0.0f)) * GameState->CubeOrientation;
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Cube);
#else

    plane Planes[4] = 
    {
        {vec2(0.0f, 1.0f), -0.4f},
        {vec2(-1.0f, 0.0f), -0.8f},
        {vec2(0.0f, -1.0f), -0.4f},
        {vec2(1.0f, 0.0f), -0.8f},
    };

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < 1;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;
        rigid_body *RigidBody = &GameObject->RigidBody;

        RigidBody->dP += dt*((RigidBody->ForceAccumulated*(1.0f / RigidBody->Mass)) + (GravityEnabled ? vec2(0.0f, -9.8f) : vec2(0.0f, 0.0f)));
        RigidBody->AngularSpeed += dt*(RigidBody->TorqueAccumulated / RigidBody->MomentOfInertia);
        RigidBody->Orientation += dt*Degrees(RigidBody->AngularSpeed);
        vec2 BodyOrientationAxisX = vec2(Cos(Radians(RigidBody->Orientation)), Sin(Radians(RigidBody->Orientation)));
        vec2 BodyOrientationAxisY = Perp(BodyOrientationAxisX);
        r32 BodyHalfWidth = 0.5f*GameObject->Width;
        r32 BodyHalfHeight = 0.5f*GameObject->Height;

        vec2 GameObjectDeltaP = dt*RigidBody->dP;
        r32 dtRemaining = dt;
        for(u32 Iteration = 0;
            Iteration < 100;
            Iteration++)
        {
            r32 GameObjectDeltaPLength = Length(GameObjectDeltaP);
            if((GameObjectDeltaPLength > 0.0f) && (dtRemaining > 0.0f))
            {
                r32 t = 1.0f;
                vec2 PointOfContact;
                vec2 CollisionNormal;
                b32 Collision = false;
                vec2 DesiredP = RigidBody->P + GameObjectDeltaP;

                for(u32 PlaneIndex = 0;
                    PlaneIndex < ArrayCount(Planes);
                    PlaneIndex++)
                {
                    plane *Plane = Planes + PlaneIndex;
                    b32 CollisionWithPlane = false;

                    r32 Radius = (BodyHalfWidth*Absolute(Dot(BodyOrientationAxisX, Plane->N)) + BodyHalfHeight*Absolute(Dot(BodyOrientationAxisY, Plane->N)));
                    r32 DistFromHeroCenterToPlane = Dot(Plane->N, RigidBody->P) - Plane->D;
                    if(Absolute(DistFromHeroCenterToPlane) < Radius)
                    {
                        // NOTE(georgy): Penetration
                        RigidBody->P += (Radius - Absolute(DistFromHeroCenterToPlane) + 0.001f)*Plane->N;
                        DistFromHeroCenterToPlane = Dot(Plane->N, RigidBody->P) - Plane->D;

                        Assert(Absolute(DistFromHeroCenterToPlane) >= Radius);
                    }   

                    {
                        r32 Denom = Dot(Plane->N, GameObjectDeltaP);
                        if((Denom * DistFromHeroCenterToPlane) >= 0.0f)
                        {
                            // NOTE(georgy): Moving parallel to or away from the plane
                            CollisionWithPlane = false;
                        }
                        else
                        {
                            r32 PlaneDisplace = (DistFromHeroCenterToPlane > 0.0f) ? Radius : -Radius;
                            r32 NewT = (PlaneDisplace - DistFromHeroCenterToPlane) / Denom;
                            if((NewT <= t) && (NewT >= 0.0f))
                            {
                                t = NewT - 0.1f;
                                if(t < 0.0f) t = NewT;
                                CollisionWithPlane = true;

                                PointOfContact = GetPointOfContact(RigidBody, t, GameObjectDeltaP, 
                                                                   BodyHalfWidth, BodyOrientationAxisX, BodyHalfHeight, BodyOrientationAxisY,
                                                                   *Plane, Radius);
                                CollisionNormal = Plane->N;
                            }
                        }
                    }

                    if(CollisionWithPlane) 
                    {
                        Collision = true;
                    }
                }

                for(u32 TestGameObjectIndex = 0;
                    TestGameObjectIndex < GameState->GameObjectCount;
                    TestGameObjectIndex++)
                {
                    if(TestGameObjectIndex != GameObjectIndex)
                    {
                        game_object *TestGameObject = GameState->GameObjects + TestGameObjectIndex;
                        
                        r32 NewT = 10000.0f;
                        if(TestIntersection(GameObject, GameObjectDeltaP, TestGameObject, &NewT))
                        {
                            if((NewT <= t) && (NewT >= 0.0f))
                            {
                                t = NewT - 0.1f;
                                PointOfContact = ????;
                                CollisionNormal = ????;
                                Collision = true;
                            }
                        }
                    }
                }

                dtRemaining -= t*dtRemaining;
                RigidBody->P += t*GameObjectDeltaP;
                GameObjectDeltaP = DesiredP - RigidBody->P;
                if(Collision)
                {
                    r32 CoeffOfRestitution = 0.25f;
                    r32 OneOverMass = 1.0f / RigidBody->Mass;
                    r32 OneOverMomentOfInertia = 1.0f / RigidBody->MomentOfInertia;

                    vec2 PointRotationDir = Perp(PointOfContact - RigidBody->P);
                    
                    vec2 Velocity = RigidBody->dP + RigidBody->AngularSpeed*PointRotationDir;

                    r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(Velocity, CollisionNormal);
                    r32 ImpulseDenom = Dot(CollisionNormal, CollisionNormal)*OneOverMass + 
                                            OneOverMomentOfInertia*Square(Dot(PointRotationDir, CollisionNormal));

                    r32 ImpulseMagnitude = ImpulseNom / ImpulseDenom;
                
                    RigidBody->dP += (ImpulseMagnitude * OneOverMass)*CollisionNormal;
                    RigidBody->AngularSpeed += Dot(PointRotationDir, ImpulseMagnitude*CollisionNormal) * OneOverMomentOfInertia;

                    GameObjectDeltaP = dtRemaining*RigidBody->dP;
                }
            }
            else
            {
                break;
            }
        }
    }

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;
        mat4 Model = Translation(vec3(GameObject->RigidBody.P, 0.0f)) * 
                     Rotation(GameObject->RigidBody.Orientation, vec3(0.0f, 0.0f, 1.0f)) * 
                     Scaling(vec3(GameObject->Width, GameObject->Height, 1.0f));
        PushMat4(RenderCommandBuffer, "Model", &Model);
        DrawModel(RenderCommandBuffer, GameObject->Model);
    }

    mat4 Model = Identity();
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawLine(RenderCommandBuffer, vec3(-1.0f, Planes[0].D*Planes[0].N.y, 0.0f), vec3(1.0f, Planes[0].D*Planes[0].N.y, 0.0f));
    DrawLine(RenderCommandBuffer, vec3(Planes[1].D*Planes[1].N.x, -1.0f, 0.0f), vec3(Planes[1].D*Planes[1].N.x, 1.0f, 0.0f));
    DrawLine(RenderCommandBuffer, vec3(-1.0f, Planes[2].D*Planes[2].N.y, 0.0f), vec3(1.0f, Planes[2].D*Planes[2].N.y, 0.0f));
    DrawLine(RenderCommandBuffer, vec3(Planes[3].D*Planes[3].N.x, -1.0f, 0.0f), vec3(Planes[3].D*Planes[3].N.x, 1.0f, 0.0f));

#endif
}
