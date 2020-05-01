#pragma once

#define PI 3.14159265358979323846f

#define Epsilon (1.19e-7f)
#define I32_MIN (-2147483647 - 1)
#define I32_MAX 2147483647
#define U32_MAX 0xFFFFFFFF
#define FLT_MAX 3.402823466e+38F

#define Max(A, B) (((A) > (B)) ? (A) : (B))
#define Min(A, B) (((A) < (B)) ? (A) : (B))

union vec2
{
    struct
    {
        r32 x, y;
    };

    struct
    {
        r32 u, v;
    };

    r32 E[2];

    vec2() {}
    vec2(r32 X, r32 Y) { x = X; y = Y; }
};

union vec3
{
    struct
    {
        r32 x, y, z;
    };

    struct
    {
        r32 u, v, w;
    };

    struct
    {
        r32 r, g, b;
    };

    struct 
    {
        vec2 xy;
        r32 Ignored_0;
    };
    struct 
    {
        r32 Ignored_1;
        vec2 yz;
    };
    struct 
    {
        vec2 uv;
        r32 Ignored_2;
    };
    struct 
    {
        r32 Ignored_3;
        vec2 vw;
    };

    r32 E[3];

    vec3() {}
    vec3(r32 X, r32 Y, r32 Z) { x = X; y = Y; z = Z; }
    vec3(vec2 XY, r32 Z) { x = XY.x; y = XY.y; z = Z; }
};

union vec4
{
    struct
    {
        r32 x, y, z, w;
    };

    struct
    {
        r32 r, g, b, a;
    };

    struct 
    {
        vec3 xyz;
        r32 Ignored_0;
    };
    struct
    {
        vec2 xy;
        r32 Ignored_1;
        r32 Ignored_2;
    };
    struct
    {
        r32 Ignored_3;
        vec2 yz;
        r32 Ignored_4;
    };
    struct
    {
        r32 Ignored_5;
        r32 Ignored_6;
        vec2 zw;
    };

    struct
    {
        vec3 rgb;
        r32 Ignored_7;
    };

    r32 E[4];

    vec4() {}
    vec4(r32 X, r32 Y, r32 Z, r32 W) { x = X; y = Y; z = Z; w = W; }
};

inline vec2
vec2i(i32 X, i32 Y)
{
    vec2 Result = vec2((r32)X, (r32)Y);

    return(Result);
}

inline vec2
vec2i(u32 X, u32 Y)
{
    vec2 Result = vec2((r32)X, (r32)Y);

    return(Result);
}

inline vec3
vec3i(i32 X, i32 Y, i32 Z)
{
    vec3 Result = vec3((r32)X, (r32)Y, (r32)Z);

    return(Result);
}

inline vec3
vec3i(u32 X, u32 Y, u32 Z)
{
    vec3 Result = vec3((r32)X, (r32)Y, (r32)Z);

    return(Result);
}

inline vec4
vec4i(i32 X, i32 Y, i32 Z, i32 W)
{
    vec4 Result = vec4((r32)X, (r32)Y, (r32)Z, (r32)W);

    return(Result);
}

inline vec4
vec4i(u32 X, u32 Y, u32 Z, u32 W)
{
    vec4 Result = vec4((r32)X, (r32)Y, (r32)Z, (r32)W);

    return(Result);
}

union mat3
{
    r32 E[9];
    struct 
    {
        r32 a11, a21, a31;
        r32 a12, a22, a32;
        r32 a13, a23, a33;
    };
};

union mat4
{
    r32 E[16];

    struct 
    {
        r32 a11, a21, a31, a41;
        r32 a12, a22, a32, a42;
        r32 a13, a23, a33, a43;
        r32 a14, a24, a34, a44;
    };
};
