#pragma once

// TODO(georgy): Get rid of math.h
#include <math.h>

inline r32
SquareRoot(r32 Value)
{
    r32 Result = sqrtf(Value);

    return(Result);
}

inline r32
Sin(r32 Angle)
{
    r32 Result = sinf(Angle);

    return(Result);
}

inline r32
Cos(r32 Angle)
{
    r32 Result = cosf(Angle);

    return(Result);
}

inline r32
Tan(r32 Angle)
{
    r32 Result = tanf(Angle);

    return(Result);
}

inline r32
ATan2(r32 Y, r32 X)
{
    r32 Result = 0.0f;

    if((Y != 0.0f) || (X != 0.0f))
    {
        // NOTE(georgy): Our Z axis (Y, here) points south when positive. 
        // But atan2f function wants opposite of that. So, we negate it here.
        Y = -Y;
        Result = atan2f(Y, X);
    }

    return(Result);
}