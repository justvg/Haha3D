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

#define Assert(Expression) if(!(Expression)) { *(int *)0 = 0; }
#define ArrayCount(Array) (sizeof(Array)/sizeof((Array)[0]))

#include "haha3d_math.h"
#include "haha3d_render_command_buffer.h"

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
};

inline b32
WasDown(button *Button)
{
    b32 Result = (Button->EndedDown && (Button->HalfTransitionCount == 1)) ||
                 (Button->HalfTransitionCount > 1);

    return(Result);
}

#define PLATFORM_INIT_BUFFERS(name) void name(void **Handle, u32 Size, r32 *Vertices)
typedef PLATFORM_INIT_BUFFERS(platform_init_buffers);

#define PLATFORM_COMPILE_SHADER(name) void name(void **Handle, char *ShaderFilepath)
typedef PLATFORM_COMPILE_SHADER(platform_compile_shader);

struct platform_api
{
    platform_init_buffers *InitBuffers;
    platform_compile_shader *CompileShader;
};

#define GAME_UPDATE_AND_RENDER(name) void name(game_input *Input, render_command_buffer *RenderCommandBuffer, u32 WindowWidth, u32 WindowHeight, platform_api PlatformAPI)
typedef GAME_UPDATE_AND_RENDER(game_update_and_render);