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
