#include "haha3d_math.h"

inline r32
Radians(r32 Angle)
{
    r32 Result = (Angle / 180.0f) * PI;;

    return(Result);
}

inline r32
Degrees(r32 Rad)
{
    r32 Result = (Rad / PI) * 180.0f;

    return(Result);
}

inline r32
Lerp(r32 A, r32 B, r32 t)
{
    r32 Result = A + (B - A)*t;

    return(Result);
}

// 
// NOTE(georgy): vec2
// 

inline vec2
operator+ (vec2 A, vec2 B)
{
    vec2 Result;

    Result.x = A.x + B.x;
    Result.y = A.y + B.y;

    return(Result);
}

inline vec2
operator- (vec2 A, vec2 B)
{
    vec2 Result;

    Result.x = A.x - B.x;
    Result.y = A.y - B.y;

    return(Result);
}

inline vec2
Hadamard(vec2 A, vec2 B)
{
    vec2 Result;

    Result.x = A.x * B.x;
    Result.y = A.y * B.y;

    return(Result);
}

inline vec2
operator* (vec2 A, r32 B)
{
    vec2 Result;

    Result.x = A.x * B;
    Result.y = A.y * B;

    return(Result);
}

inline vec2
operator* (r32 B, vec2 A)
{
    vec2 Result = A * B;

    return(Result);
}

inline vec2 &
operator+= (vec2 &A, vec2 B)
{
    A = A + B;
    
    return(A);
}

inline vec2 &
operator-= (vec2 &A, vec2 B)
{
    A = A - B;
    
    return(A);
}

inline vec2 &
operator*= (vec2 &A, r32 B)
{
    A = A * B;
    
    return(A);
}

inline vec2
operator- (vec2 A)
{
    vec2 Result;

    Result.x = -A.x;
    Result.y = -A.y;
    
    return(Result);
}

inline r32
Dot(vec2 A, vec2 B)
{
    r32 Result = A.x*B.x + A.y*B.y;

    return(Result);
}

inline r32
LengthSq(vec2 A)
{
    r32 Result = Dot(A, A);

    return(Result);
}

inline r32
Length(vec2 A)
{
    r32 Result = SquareRoot(LengthSq(A));

    return(Result);
}

inline vec2
Normalize(vec2 A)
{
    vec2 Result = A * (1.0f / Length(A));

    return(Result);
}

inline vec2
Perp(vec2 A)
{
    vec2 Result = vec2(-A.y, A.x);
    
    return(Result);
}

inline r32
Cross2D(vec2 A, vec2 B)
{
    r32 Result = Dot(Perp(A), B);

    return(Result);
}

inline vec2
Lerp(vec2 A, vec2 B, r32 t)
{
    vec2 Result = A + (B - A)*t;

    return(Result);
}

// 
// NOTE(georgy): vec3
// 

inline vec3
operator+ (vec3 A, vec3 B)
{
    vec3 Result;

    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    Result.z = A.z + B.z;

    return(Result);
}

inline vec3
operator- (vec3 A, vec3 B)
{
    vec3 Result;

    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;

    return(Result);
}

inline vec3
Hadamard(vec3 A, vec3 B)
{
    vec3 Result;

    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;

    return(Result);
}

inline vec3
operator* (vec3 A, r32 B)
{
    vec3 Result;

    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;

    return(Result);
}

inline vec3
operator* (r32 B, vec3 A)
{
    vec3 Result = A * B;

    return(Result);
}

inline vec3 &
operator+= (vec3 &A, vec3 B)
{
    A = A + B;
    
    return(A);
}

inline vec3 &
operator-= (vec3 &A, vec3 B)
{
    A = A - B;
    
    return(A);
}

inline vec3 &
operator*= (vec3 &A, r32 B)
{
    A = A * B;
    
    return(A);
}

inline vec3
operator- (vec3 A)
{
    vec3 Result;

    Result.x = -A.x;
    Result.y = -A.y;
    Result.z = -A.z;
    
    return(Result);
}

inline r32
Dot(vec3 A, vec3 B)
{
    r32 Result = A.x*B.x + A.y*B.y + A.z*B.z;

    return(Result);
}

inline r32
LengthSq(vec3 A)
{
    r32 Result = Dot(A, A);

    return(Result);
}

inline r32
Length(vec3 A)
{
    r32 Result = SquareRoot(LengthSq(A));

    return(Result);
}

inline vec3
Normalize(vec3 A)
{
    vec3 Result = A * (1.0f / Length(A));

    return(Result);
}

inline vec3
Cross(vec3 A, vec3 B)
{
    vec3 Result;

    Result.x = A.y*B.z - B.y*A.z;
    Result.y = A.z*B.x - B.z*A.x;
    Result.z = A.x*B.y - B.x*A.y;

    return(Result);    
}

inline vec3
Lerp(vec3 A, vec3 B, r32 t)
{
    vec3 Result = A + (B - A)*t;

    return(Result);
}


// 
// NOTE(georgy): vec4
// 

inline vec4
operator+ (vec4 A, vec4 B)
{
    vec4 Result;

    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    Result.z = A.z + B.z;
    Result.w = A.w + B.w;

    return(Result);
}

inline vec4
operator- (vec4 A, vec4 B)
{
    vec4 Result;

    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    Result.w = A.w - B.w;

    return(Result);
}

inline vec4
Hadamard(vec4 A, vec4 B)
{
    vec4 Result;

    Result.x = A.x * B.x;
    Result.y = A.y * B.y;
    Result.z = A.z * B.z;
    Result.w = A.w * B.w;

    return(Result);
}

inline vec4
operator* (vec4 A, r32 B)
{
    vec4 Result;

    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    Result.w = A.w * B;

    return(Result);
}

inline vec4
operator* (r32 B, vec4 A)
{
    vec4 Result = A * B;

    return(Result);
}

inline vec4 &
operator+= (vec4 &A, vec4 B)
{
    A = A + B;
    
    return(A);
}

inline vec4 &
operator-= (vec4 &A, vec4 B)
{
    A = A - B;
    
    return(A);
}

inline vec4 &
operator*= (vec4 &A, r32 B)
{
    A = A * B;
    
    return(A);
}

inline vec4
operator- (vec4 A)
{
    vec4 Result;

    Result.x = -A.x;
    Result.y = -A.y;
    Result.z = -A.z;
    Result.w = -A.w;
    
    return(Result);
}

inline r32
Dot(vec4 A, vec4 B)
{
    r32 Result = A.x*B.x + A.y*B.y + A.z*B.z + A.w*B.w;

    return(Result);
}

inline r32
LengthSq(vec4 A)
{
    r32 Result = Dot(A, A);

    return(Result);
}

inline r32
Length(vec4 A)
{
    r32 Result = SquareRoot(LengthSq(A));

    return(Result);
}

inline vec4
Normalize(vec4 A)
{
    vec4 Result = A * (1.0f / Length(A));

    return(Result);
}

inline vec4
Lerp(vec4 A, vec4 B, r32 t)
{
    vec4 Result = A + (B - A)*t;

    return(Result);
}

// 
// NOTE(georgy): plane
// 

struct plane
{
    vec2 N;
    r32 D;
};

// 
// NOTE(georgy): mat3
// 

inline mat3
Identity3x3(r32 Diagonal = 1.0f)
{
    mat3 Result = {};

    Result.a11 = Diagonal;
    Result.a22 = Diagonal;
    Result.a33 = Diagonal;

    return(Result);
}

inline mat3
Scaling3x3(r32 Scale)
{
    mat3 Result = {};

    Result.a11 = Scale;
    Result.a22 = Scale;
    Result.a33 = Scale;

    return(Result);
}

inline mat3
Scaling3x3(vec3 Scale)
{
    mat3 Result = {};

    Result.a11 = Scale.x;
    Result.a22 = Scale.y;
    Result.a33 = Scale.z;

    return(Result);
}

internal mat3
Rotation3x3(r32 Angle, vec3 Axis)
{
    mat3 Result;

    r32 Rad = Radians(Angle);
    r32 Sine = Sin(Rad);
    r32 Cosine = Cos(Rad);

    Axis = Normalize(Axis);

    Result.a11 = Axis.x*Axis.x*(1.0f - Cosine) + Cosine;
	Result.a21 = Axis.x*Axis.y*(1.0f - Cosine) + Axis.z*Sine;
	Result.a31 = Axis.x*Axis.z*(1.0f - Cosine) - Axis.y*Sine;

	Result.a12 = Axis.x*Axis.y*(1.0f - Cosine) - Axis.z*Sine;
	Result.a22 = Axis.y*Axis.y*(1.0f - Cosine) + Cosine;
	Result.a32 = Axis.y*Axis.z*(1.0f - Cosine) + Axis.x*Sine;

	Result.a13 = Axis.x*Axis.z*(1.0f - Cosine) + Axis.y*Sine;
	Result.a23 = Axis.y*Axis.z*(1.0f - Cosine) - Axis.x*Sine;
	Result.a33 = Axis.z*Axis.z*(1.0f - Cosine) + Cosine;

    return(Result);
}

internal mat3
operator*(mat3 A, mat3 B)
{
    mat3 Result;

    for(u32 Row = 0;
        Row < 3;
        Row++)
    {
        for(u32 Column = 0;
            Column < 3;
            Column++)
        {
            r32 Sum = 0.0f;
            for(u32 E = 0;
                E < 3;
                E++)
            {
                Sum += A.E[Row + E*3] * B.E[Column*3 + E];
            }
            Result.E[Row + Column*3] = Sum;
        }
    }

    return(Result);
}

// 
// NOTE(georgy): mat4
// 

inline mat4
Identity(r32 Diagonal = 1.0f)
{
    mat4 Result = {};

    Result.a11 = Diagonal;
    Result.a22 = Diagonal;
    Result.a33 = Diagonal;
    Result.a44 = Diagonal;

    return(Result);
}

inline mat4
Translation(vec3 Translate)
{
    mat4 Result = Identity(1.0f);

    Result.a14 = Translate.x;
    Result.a24 = Translate.y;
    Result.a34 = Translate.z;

    return(Result);
}

inline mat4
Scaling(r32 Scale)
{
    mat4 Result = {};

    Result.a11 = Scale;
    Result.a22 = Scale;
    Result.a33 = Scale;
    Result.a44 = 1.0f;

    return(Result);
}

inline mat4
Scaling(vec3 Scale)
{
    mat4 Result = {};

    Result.a11 = Scale.x;
    Result.a22 = Scale.y;
    Result.a33 = Scale.z;
    Result.a44 = 1.0f;

    return(Result);
}

internal mat4
Rotation(r32 Angle, vec3 Axis)
{
    mat4 Result;

    r32 Rad = Radians(Angle);
    r32 Sine = Sin(Rad);
    r32 Cosine = Cos(Rad);

    Axis = Normalize(Axis);

    Result.a11 = Axis.x*Axis.x*(1.0f - Cosine) + Cosine;
	Result.a21 = Axis.x*Axis.y*(1.0f - Cosine) + Axis.z*Sine;
	Result.a31 = Axis.x*Axis.z*(1.0f - Cosine) - Axis.y*Sine;
	Result.a41 = 0.0f;

    Result.a12 = Axis.x*Axis.y*(1.0f - Cosine) - Axis.z*Sine;
	Result.a22 = Axis.y*Axis.y*(1.0f - Cosine) + Cosine;
	Result.a32 = Axis.y*Axis.z*(1.0f - Cosine) + Axis.x*Sine;
	Result.a42 = 0.0f;

    Result.a13 = Axis.x*Axis.z*(1.0f - Cosine) + Axis.y*Sine;
	Result.a23 = Axis.y*Axis.z*(1.0f - Cosine) - Axis.x*Sine;
	Result.a33 = Axis.z*Axis.z*(1.0f - Cosine) + Cosine;
	Result.a43 = 0.0f;

	Result.a14 = 0.0f;
	Result.a24 = 0.0f;
	Result.a34 = 0.0f;
	Result.a44 = 1.0f;

    return(Result);
}

internal mat4
LookAt(vec3 From, vec3 Target, vec3 UpAxis = vec3(0.0f, 1.0f, 0.0f))
{
    vec3 Forward = Normalize(From - Target);
    vec3 Right = Normalize(Cross(UpAxis, Forward));
    vec3 Up = Cross(Forward, Right);

    mat4 Result;

    Result.a11 = Right.x;
    Result.a21 = Up.x;
    Result.a31 = Forward.x;
    Result.a41 = 0.0f;

    Result.a12 = Right.y;
    Result.a22 = Up.y;
    Result.a32 = Forward.y;
    Result.a42 = 0.0f;

    Result.a13 = Right.z;
    Result.a23 = Up.z;
    Result.a33 = Forward.z;
    Result.a43 = 0.0f;

    Result.a14 = -Dot(Right, From);
    Result.a24 = -Dot(Up, From);
    Result.a34 = -Dot(Forward, From);
    Result.a44 = 1.0f;

    return(Result);
}

internal mat4
ViewRotationMatrixFromDirection(vec3 Dir, vec3 UpAxis = vec3(0.0f, 1.0f, 0.0f))
{
    vec3 Forward = Normalize(-Dir);
    vec3 Right = Normalize(Cross(UpAxis, Forward));
    vec3 Up = Cross(Forward, Right);

    mat4 Result;

    Result.a11 = Right.x;
    Result.a21 = Up.x;
    Result.a31 = Forward.x;
    Result.a41 = 0.0f;

    Result.a12 = Right.y;
    Result.a22 = Up.y;
    Result.a32 = Forward.y;
    Result.a42 = 0.0f;

    Result.a13 = Right.z;
    Result.a23 = Up.z;
    Result.a33 = Forward.z;
    Result.a43 = 0.0f;

    Result.a14 = 0.0f;
    Result.a24 = 0.0f;
    Result.a34 = 0.0f;
    Result.a44 = 1.0f;

    return(Result);
}

internal mat4 
Ortho(r32 Bottom, r32 Top, r32 Left, r32 Right, r32 Near, r32 Far)
{
    mat4 Result = {};

    Result.a11 = 2.0f / (Right - Left);
    Result.a22 = 2.0f / (Top - Bottom);
    Result.a33 = -2.0f / (Far - Near);
    Result.a14 = -(Right + Left) / (Right - Left);
    Result.a24 = -(Top + Bottom) / (Top - Bottom);
	Result.a34 = -(Far + Near) / (Far - Near);
	Result.a44 = 1.0f;

    return(Result);
}

internal mat4
Perspective(r32 FoV, r32 AspectRatio, r32 Near, r32 Far)
{
    r32 Scale = Tan(Radians(FoV) * 0.5f) * Near;
    r32 Top = Scale;
    r32 Bottom = -Top;
    r32 Right = Scale * AspectRatio;
    r32 Left = -Right;

    mat4 Result = {};

    Result.a11 = 2.0f * Near / (Right - Left);
	Result.a22 = 2.0f * Near / (Top - Bottom);
	Result.a13 = (Right + Left) / (Right - Left);
	Result.a23 = (Top + Bottom) / (Top - Bottom);
	Result.a33 = -(Far + Near) / (Far - Near);
	Result.a43 = -1.0f;
	Result.a34 = -(2.0f * Far*Near) / (Far - Near);

    return(Result);
}

internal mat4
operator*(mat4 A, mat4 B)
{
    mat4 Result;

    for(u32 Row = 0;
        Row < 4;
        Row++)
    {
        for(u32 Column = 0;
            Column < 4;
            Column++)
        {
            r32 Sum = 0.0f;
            for(u32 E = 0;
                E < 4;
                E++)
            {
                Sum += A.E[Row + E*4] * B.E[Column*4 + E];
            }
            Result.E[Row + Column*4] = Sum;
        }
    }

    return(Result);
}
