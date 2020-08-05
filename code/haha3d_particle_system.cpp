#include <vector>

template <typename T>
struct animated_value_key
{
    r32 Time;
    T Value;
};

template <typename T>
struct animated_value
{
    std::vector<animated_value_key<T>> Keys;

    void AddKey(r32 Time, T Value);
    T GetValueForTime(r32 Time);
};

template <typename T>
void animated_value<T>::AddKey(r32 Time, T Value)
{
    animated_value_key<T> Key = { Time, Value };
    Keys.push_back(Key);
}

template <typename T>
T animated_value<T>::GetValueForTime(r32 Time)
{
    // NOTE(georgy): We need to have at least 2 keys
    Assert(Keys.size() > 2);
    Assert((Time >= 0.0f) && (Time <= 1.0f));

    u32 FirstKey = 0;
    u32 SecondKey = 0;
    for(u32 KeyIndex = 0;
        KeyIndex < Keys.size() - 1;
        KeyIndex++)
    {
        if(Time <= Keys[KeyIndex + 1].Time)
        {
            FirstKey = KeyIndex;
            SecondKey = KeyIndex + 1;
            break;
        }
    }

    r32 TimeBetweenKeys = Keys[SecondKey].Time - Keys[FirstKey].Time;
    r32 t = (Time - Keys[FirstKey].Time) / TimeBetweenKeys;

    T Result = Lerp(Keys[FirstKey].Value, Keys[SecondKey].Value, t);

    return(Result);
}

struct particle
{
    vec3 P;
    vec3 Velocity;
    vec3 Acceleration;

    r32 Rotation;
    r32 RotationSpeed;

    vec2 ScaleVarianceStart;
    vec2 Scale;

    vec4 ColorMult;

    r32 CurrentLifeTime;
    r32 LifeTime;
};

struct particle_emmiter
{
    vec3 P;
    vec3 SpawnArea;

    vec3 InitialVelocity;
    vec3 VelocityVariance;

    vec3 InitialAcceleration;

    r32 InitialRotation;
    r32 RotationVariance;
    r32 RotationSpeedVariance;

    vec2 InitialScale;
    animated_value<vec2> Scale;
    r32 ScaleVariance;

    vec4 InitialColorMult;
    vec4 ColorMultVariance;

    r32 ParticleLifeTime;
    r32 ParticleLifeTimeVariance;

    r32 PartialParticle;
    u32 MaxParticlesCount;
    std::vector<particle> Particles;
};
