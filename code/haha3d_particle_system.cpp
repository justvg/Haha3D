#include "haha3d_particle_system.h"

internal void
UpdateParticles(particle_emmiter *Emitter, r32 dt)
{
    // NOTE(georgy): Update particles
    for(u32 ParticleIndex = 0;
        ParticleIndex < Emitter->Particles.size();
        )
    {
        particle *Particle = &Emitter->Particles[0] + ParticleIndex;

        Particle->CurrentLifeTime += dt;

        if(Particle->CurrentLifeTime > Particle->LifeTime)
        {
            Emitter->Particles[ParticleIndex] = Emitter->Particles[Emitter->Particles.size() - 1];
            Emitter->Particles.pop_back();
        }
        else
        {
            r32 t = Particle->CurrentLifeTime / Particle->LifeTime;

            Particle->AccumulatedAcceleration += dt*Emitter->InitialAcceleration;
            Particle->Velocity = Particle->InitialVelocityVariance + Particle->AccumulatedAcceleration + Emitter->Velocity.GetValueForTime(t);
            Particle->P += dt*Particle->Velocity;

            Particle->Rotation += dt*(Particle->InitialRotationSpeedVariance + Emitter->RotationSpeed.GetValueForTime(t));
            Particle->Scale = Particle->InitialScaleVariance + Emitter->Scale.GetValueForTime(t);

            Particle->ColorMult = Particle->InitialColorMultVariance + Emitter->ColorMult.GetValueForTime(t);

            ParticleIndex++;
        }
    }

    // NOTE(georgy): Add new particles
    Emitter->CurrentLifeTime += dt;
    if(Emitter->CurrentLifeTime > Emitter->LifeTime)
    {
        if(Emitter->Looped)
        {
            Emitter->CurrentLifeTime -= Emitter->LifeTime;
        }
        else
        {
            // TODO(georgy): Delete emitter here
        }
    }
    r32 t = Emitter->CurrentLifeTime / Emitter->LifeTime;
    r32 NewParticlesCountReal = Emitter->EmissionPerSecond.GetValueForTime(t)*dt;
    u32 NewParticlesCount = (u32)NewParticlesCountReal;
    
    Emitter->PartialParticle += NewParticlesCountReal - NewParticlesCount;
    if(Emitter->PartialParticle > 1.0f)
    {
        NewParticlesCount++;
        Emitter->PartialParticle -= 1.0f;
    }

    for(u32 NewParticle = 0;
        NewParticle < NewParticlesCount;
        NewParticle++)
    {
        if(Emitter->Particles.size() < Emitter->MaxParticlesCount)
        {
            particle Particle;

            vec3 Offset = Hadamard(vec3(Emitter->RandDistribution(Emitter->RandGenerator), Emitter->RandDistribution(Emitter->RandGenerator), Emitter->RandDistribution(Emitter->RandGenerator)),
                                   0.5f*Emitter->SpawnArea);
            Particle.P = Emitter->P + Offset;

            Particle.InitialVelocityVariance = Hadamard(vec3(Emitter->RandDistribution(Emitter->RandGenerator), Emitter->RandDistribution(Emitter->RandGenerator), Emitter->RandDistribution(Emitter->RandGenerator)),
                                                        Emitter->VelocityVariance);
            Particle.Velocity = Emitter->Velocity.GetValueForTime(0.0f) + Particle.InitialVelocityVariance;

            Particle.AccumulatedAcceleration = vec3(0.0f, 0.0f, 0.0f);

            Particle.Rotation = Emitter->InitialRotation + Emitter->RandDistribution(Emitter->RandGenerator)*Emitter->RotationVariance;
            Particle.InitialRotationSpeedVariance = Emitter->RandDistribution(Emitter->RandGenerator)*Emitter->RotationSpeedVariance;
            Particle.RotationSpeed = Emitter->RotationSpeed.GetValueForTime(0.0f) + Particle.InitialRotationSpeedVariance;

            Particle.InitialScaleVariance = Emitter->RandDistribution(Emitter->RandGenerator)*vec2(Emitter->ScaleVariance, Emitter->ScaleVariance); 
            Particle.Scale = Emitter->Scale.GetValueForTime(0.0f) + Particle.InitialScaleVariance;

            Particle.InitialColorMultVariance = Hadamard(vec4(Emitter->RandDistribution(Emitter->RandGenerator), Emitter->RandDistribution(Emitter->RandGenerator), Emitter->RandDistribution(Emitter->RandGenerator), Emitter->RandDistribution(Emitter->RandGenerator)),
                                                              Emitter->ColorMultVariance);
            Particle.ColorMult = Emitter->ColorMult.GetValueForTime(0.0f) + Particle.InitialColorMultVariance;

            Particle.CurrentLifeTime = 0.0f;
            Particle.LifeTime = Emitter->ParticleLifeTime + Emitter->RandDistribution(Emitter->RandGenerator)*Emitter->ParticleLifeTimeVariance;

            Emitter->Particles.push_back(Particle);
        }
    }
}

internal void
RenderParticles(particle_emmiter *Emitter, shader Shader, vec3 CameraForward, vec3 CameraRight, u32 VertexCount)
{
    for(u32 ParticleIndex = 0;
        ParticleIndex < Emitter->Particles.size();
        ParticleIndex++)
    {
        particle *Particle = &Emitter->Particles[0] + ParticleIndex;

        vec3 CamRight = Rotation3x3(Particle->Rotation, vec3(-CameraForward)) * CameraRight;
        vec3 CamUp = Cross(-CameraForward, CamRight);

        Shader.SetVec3("WorldP", Particle->P);
        Shader.SetVec3("CameraRight", CamRight);
        Shader.SetVec3("CameraUp", CamUp);
        Shader.SetVec2("Scale", Particle->Scale);
        Shader.SetVec4("ColorMult", Particle->ColorMult);
        glDrawArrays(GL_TRIANGLES, 0, VertexCount);
    }
}