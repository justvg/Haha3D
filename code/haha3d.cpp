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
    r32 ClosestDistance = Absolute(Dot(RigidBodyCorners[0], Plane.N) - Plane.D);

    ClosestIndecies[ClosestCount++] = TheClosestToPlaneIndex;

    for(u32 Corner = 1;
        Corner < 4;
        Corner++)
    {
        r32 Distance = Absolute(Dot(RigidBodyCorners[Corner], Plane.N) - Plane.D);
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

    Assert(ClosestCount <= 2);

    vec2 PointOfContact;
    if(ClosestCount == 2)
    {
        PointOfContact = Lerp(RigidBodyCorners[ClosestIndecies[0]], RigidBodyCorners[ClosestIndecies[1]], 0.5f);
    }
    else
    {
        PointOfContact = RigidBodyCorners[TheClosestToPlaneIndex];
    }

    return(PointOfContact);
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

        GameState->HeroQuad.RigidBody.Mass = 3.0f;
        GameState->HeroQuad.RigidBody.MomentOfInertia = (1.0f / 12.0f) * GameState->HeroQuad.RigidBody.Mass * (0.05f*0.05f+ 0.1f*0.1f);

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

    GameState->HeroQuad.RigidBody.ForceAccumulated = vec2(0.0f, 0.0f);
    GameState->HeroQuad.RigidBody.TorqueAccumulated = 0.0f;

    GameState->HeroQuad.RigidBody.ForceAccumulated += (1.0f / GameState->HeroQuad.RigidBody.Mass) * vec2(0.0f, -9.8f);

    static b32 ForceActsOnCM = true;
    if(WasDown(&Input->MouseLeft))
    {
        ForceActsOnCM = !ForceActsOnCM;
    }

    rigid_body *HeroRigidBody = &GameState->HeroQuad.RigidBody;
    r32 HeroWidth = 0.1f;
    r32 HeroHeight = 0.05f;
    r32 HeroHalfWidth = 0.5f*HeroWidth;
    r32 HeroHalfHeight = 0.5f*HeroHeight;
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

    plane Plane = {vec2(0.0f, 1.0f), -0.4f};

    HeroRigidBody->dP += dt*((HeroRigidBody->ForceAccumulated*(1.0f / HeroRigidBody->Mass)));
    HeroRigidBody->AngularSpeed += dt*(HeroRigidBody->TorqueAccumulated / HeroRigidBody->MomentOfInertia);
    HeroRigidBody->Orientation += dt*Degrees(HeroRigidBody->AngularSpeed);
    HeroOrientationAxisX = vec2(Cos(Radians(HeroRigidBody->Orientation)), Sin(Radians(HeroRigidBody->Orientation)));
    HeroOrientationAxisY = Perp(HeroOrientationAxisX);
    r32 Radius = (HeroHalfWidth*Absolute(Dot(HeroOrientationAxisX, Plane.N)) + HeroHalfHeight*Absolute(Dot(HeroOrientationAxisY, Plane.N)));

    vec2 HeroDeltaP = dt*HeroRigidBody->dP;
    r32 dtRemaining = dt;
    for(u32 Iteration = 0;
        Iteration < 100;
        Iteration++)
    {
        r32 HeroDeltaPLength = Length(HeroDeltaP);
        if((HeroDeltaPLength > 0.0f) && (dtRemaining > 0.0f))
        {
            r32 t = 1.0f;
            vec2 PointOfContact;
            b32 Collision = false;
            vec2 DesiredP = HeroRigidBody->P + HeroDeltaP;

            r32 DistFromHeroCenterToPlane = Dot(Plane.N, HeroRigidBody->P) - Plane.D;
            if(Absolute(DistFromHeroCenterToPlane) < Radius)
            {
                // NOTE(georgy): Penetration
                HeroRigidBody->P += (Radius - Absolute(DistFromHeroCenterToPlane) + 0.001f)*Plane.N;
                DistFromHeroCenterToPlane = Dot(Plane.N, HeroRigidBody->P) - Plane.D;

                Assert(Absolute(DistFromHeroCenterToPlane) >= Radius);
            }   

            {
                r32 Denom = Dot(Plane.N, HeroDeltaP);
                if((Denom * DistFromHeroCenterToPlane) >= 0.0f)
                {
                    // NOTE(georgy): Moving parallel to or away from the plane
                    Collision = false;
                }
                else
                {
                    r32 PlaneDisplace = (DistFromHeroCenterToPlane > 0.0f) ? Radius : -Radius;
                    r32 NewT = (PlaneDisplace - DistFromHeroCenterToPlane) / Denom;
                    if((NewT <= t) && (NewT >= 0.0f))
                    {
                        t = NewT - 0.1f;
                        if(t < 0.0f) t = NewT;
                        Collision = true;

                        PointOfContact = GetPointOfContact(HeroRigidBody, t, HeroDeltaP, 
                                                        HeroHalfWidth, HeroOrientationAxisX, HeroHalfHeight, HeroOrientationAxisY,
                                                        Plane, Radius);
                    }
                }
            }

            dtRemaining -= t*dtRemaining;
            HeroRigidBody->P += t*HeroDeltaP;
            HeroDeltaP = DesiredP - HeroRigidBody->P;
            if(Collision)
            {
                r32 CoeffOfRestitution = 0.25f;
                r32 OneOverMass = 1.0f / HeroRigidBody->Mass;
                r32 OneOverMomentOfInertia = 1.0f / HeroRigidBody->MomentOfInertia;

                vec2 PointRotationDir = Perp(PointOfContact - HeroRigidBody->P);
                
                vec2 Velocity = HeroRigidBody->dP + HeroRigidBody->AngularSpeed*PointRotationDir;

                r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(Velocity, Plane.N);
                r32 ImpulseDenom = Dot(Plane.N, Plane.N)*OneOverMass + 
                                        OneOverMomentOfInertia*Square(Dot(PointRotationDir, Plane.N));

                r32 ImpulseMagnitude = ImpulseNom / ImpulseDenom;
            
                HeroRigidBody->dP += (ImpulseMagnitude * OneOverMass)*Plane.N;
                HeroRigidBody->AngularSpeed += Dot(PointRotationDir, ImpulseMagnitude*Plane.N) * OneOverMomentOfInertia;

                HeroDeltaP = dtRemaining*HeroRigidBody->dP;
            }
        }
    }

    mat4 Model = Translation(vec3(HeroRigidBody->P, 0.0f)) * Rotation(HeroRigidBody->Orientation, vec3(0.0f, 0.0f, 1.0f)) * Scaling(vec3(HeroWidth, HeroHeight, 1.0f));
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Quad);

    Model = Identity();
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawLine(RenderCommandBuffer, vec3(-1.0f, Plane.D, 0.0f), vec3(1.0f, -0.4f, 0.0f));

#endif
}
