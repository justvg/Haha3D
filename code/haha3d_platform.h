#pragma once

#define global_variable static
#define internal static

#include <stdint.h>

typedef  uint8_t u8;
typedef  uint8_t b8;
typedef  uint16_t u16;
typedef  uint32_t u32;
typedef  uint64_t u64;

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef int32_t b32;

typedef float r32;
typedef double r64;

#define Kilobytes(Value) ((Value) * 1024LL)
#define Megabytes(Value) (Kilobytes(Value) * 1024LL)
#define Gigabytes(Value) (Gigabytes(Value) * 1024LL)

#define Assert(Expression) if(!(Expression)) { *(int *)0 = 0; }
#define ArrayCount(Array) (sizeof(Array)/sizeof((Array)[0]))

struct button
{
    b32 EndedDown;
    u32 HalfTransitionCount;
};

struct game_input
{
    i32 MouseX, MouseY;
    i32 MouseXDisplacement, MouseYDisplacement;

    union
    {
        button Buttons[8];
        struct
        {
            button MoveForward;
            button MoveBack;
            button MoveRight;
            button MoveLeft;

            button MouseLeft, MouseRight;

            button F4;
            button Alt;
        };
    };

    r32 dtForFrame;
};

inline b32
WasDown(button *Button)
{
    b32 Result = (Button->EndedDown && (Button->HalfTransitionCount == 1)) ||
                 (Button->HalfTransitionCount > 1);

    return(Result);
}
struct game_memory
{
    u64 PermanentStorageSize;
    void *PermanentStorage;
};