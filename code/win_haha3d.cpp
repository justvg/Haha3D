#include <Windows.h>
#include <stdint.h>

#include "glew\glew.h"
#include "glew\wglew.h"

#include <string>
#include <fstream>

#define global_variable static
#define internal static

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

#include "haha3d_intrinsics.h"
#include "haha3d_math.h"
#include "haha3d_renderer_opengl.cpp"

global_variable b32 GlobalRunning;

internal void
WinSetPixelFormat(HDC WindowDC)
{
    int SuggestedPixelFormatIndex = 0;

    // NOTE(georgy): wglChoosePixelFormatARB is initialized by glewInit() call
    if(wglChoosePixelFormatARB)
    {
        int PixelFormatAttributes[] = 
        {
            WGL_DRAW_TO_WINDOW_ARB, GL_TRUE,
            WGL_ACCELERATION_ARB, WGL_FULL_ACCELERATION_ARB,
            WGL_SUPPORT_OPENGL_ARB, GL_TRUE,
            WGL_DOUBLE_BUFFER_ARB, GL_TRUE,
            WGL_PIXEL_TYPE_ARB, WGL_TYPE_RGBA_ARB,
            WGL_COLOR_BITS_ARB, 32,
            WGL_DEPTH_BITS_ARB, 24,
            WGL_STENCIL_BITS_ARB, 8,
            WGL_SAMPLE_BUFFERS_ARB, 1,
            WGL_SAMPLES_ARB, 1,
            WGL_FRAMEBUFFER_SRGB_CAPABLE_ARB, GL_TRUE,
            0
        };

        u32 ExtendedChoose;
        wglChoosePixelFormatARB(WindowDC, PixelFormatAttributes, 0, 1, &SuggestedPixelFormatIndex, &ExtendedChoose);
    }
    else
    {
        PIXELFORMATDESCRIPTOR DesiredPixelFormat = {};
        DesiredPixelFormat.nSize = sizeof(DesiredPixelFormat);
        DesiredPixelFormat.nVersion = 1;
        DesiredPixelFormat.iPixelType = PFD_TYPE_RGBA;
        DesiredPixelFormat.dwFlags = PFD_SUPPORT_OPENGL | PFD_DRAW_TO_WINDOW | PFD_DOUBLEBUFFER;
        DesiredPixelFormat.cColorBits = 32;
        DesiredPixelFormat.cAlphaBits = 8;
        DesiredPixelFormat.cDepthBits = 24;
        DesiredPixelFormat.cStencilBits = 8;
        DesiredPixelFormat.iLayerType = PFD_MAIN_PLANE;

        SuggestedPixelFormatIndex = ChoosePixelFormat(WindowDC, &DesiredPixelFormat);
    }

    PIXELFORMATDESCRIPTOR SuggestedPixelFormat;
    DescribePixelFormat(WindowDC, SuggestedPixelFormatIndex, sizeof(SuggestedPixelFormat), &SuggestedPixelFormat);
    SetPixelFormat(WindowDC, SuggestedPixelFormatIndex, &SuggestedPixelFormat);
}

internal void
WinLoadWGLExtensionsAndInitGLEW(void)
{
    WNDCLASS WindowClass = {};
    WindowClass.lpfnWndProc = DefWindowProc;
    WindowClass.hInstance = GetModuleHandle(0);
    WindowClass.lpszClassName = "FakeWindowClass";
    
    if(RegisterClass(&WindowClass))
    {
        HWND Window = CreateWindowEx(0, WindowClass.lpszClassName, "FakeWindow",
                                     WS_OVERLAPPEDWINDOW, 
                                     CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
                                     0, 0, WindowClass.hInstance, 0);

        HDC WindowDC = GetDC(Window);
        WinSetPixelFormat(WindowDC);

        HGLRC OpenGLRC = wglCreateContext(WindowDC);
        if(wglMakeCurrent(WindowDC, OpenGLRC))
        {
            glewInit();
            wglMakeCurrent(0, 0);
        }

        wglDeleteContext(OpenGLRC);
        ReleaseDC(Window, WindowDC);
        DestroyWindow(Window);
    }
}

internal b32
WinInitOpenGL(HWND Window)
{
    b32 ModernOpenGLInitialized = false;

    WinLoadWGLExtensionsAndInitGLEW();

    HDC WindowDC = GetDC(Window);
    WinSetPixelFormat(WindowDC);

    HGLRC OpenGLRC = 0;    
    if(wglCreateContextAttribsARB)
    {
        int ContextAttributes[] = 
        {
            WGL_CONTEXT_MAJOR_VERSION_ARB, 3,
            WGL_CONTEXT_MINOR_VERSION_ARB, 3,
            0
        };

        OpenGLRC = wglCreateContextAttribsARB(WindowDC, 0, ContextAttributes);
    }

    if(OpenGLRC)
    {
        if(wglMakeCurrent(WindowDC, OpenGLRC))
        {
            ModernOpenGLInitialized = true;

            InitOpenGLProperties();
        }
    }

    ReleaseDC(Window, WindowDC);

    return(ModernOpenGLInitialized);
}

LRESULT CALLBACK
WinWindowCallback(HWND Window, UINT Message, WPARAM WParam, LPARAM LParam)
{
    LRESULT Result = 0;

    switch(Message)
    {
        case WM_ACTIVATEAPP:
        {
            if(WParam)
            {
                ShowCursor(FALSE);
            }
            else
            {
                ShowCursor(TRUE);
            }
        } break;

        case WM_CLOSE:
        case WM_DESTROY:
        {
            GlobalRunning = false;
        } break;

		case WM_PAINT:
		{
			PAINTSTRUCT Paint;
			HDC DeviceContext = BeginPaint(Window, &Paint);

			EndPaint(Window, &Paint);
		} break;

        default:
        {
            Result = DefWindowProc(Window, Message, WParam, LParam);
        }
    }

    return(Result);
}

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

inline void
WinProcessKey(button *Button, b32 IsDown)
{
    if(Button->EndedDown != IsDown)
    {  
        Button->HalfTransitionCount++;
    }

    Button->EndedDown = IsDown;
}

int CALLBACK
WinMain(HINSTANCE Instance, HINSTANCE PrevInstance, LPSTR CommandLine, int ShowCode)
{
    WNDCLASS WindowClass = {};
    WindowClass.style = CS_VREDRAW | CS_HREDRAW;
    WindowClass.lpfnWndProc = WinWindowCallback;
    WindowClass.hInstance = Instance;
    WindowClass.hCursor = LoadCursor(0, IDC_ARROW);
    WindowClass.lpszClassName = "Haha3DWindowClass";
    WindowClass.hbrBackground = CreateSolidBrush(RGB(255, 0, 0));

    if(RegisterClass(&WindowClass))
    {
        HWND Window = CreateWindowEx(0, WindowClass.lpszClassName,
                                     "Haha3D", WS_OVERLAPPEDWINDOW,
                                     CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
                                     0, 0, WindowClass.hInstance, 0);

        if(Window)
        {
            GlobalRunning = true;

            if(!WinInitOpenGL(Window))
            {
                GlobalRunning = false;
                MessageBox(Window, "Can't initialize OpenGL.", 0, MB_OK);
            }
            else
            {
                ShowWindow(Window, SW_SHOW);
            }
            
            RAWINPUTDEVICE RIDs[2];
            RIDs[0].usUsagePage = 0x01;
            RIDs[0].usUsage = 0x06; // NOTE(georgy): keyboard input
            RIDs[0].dwFlags = RIDEV_NOLEGACY;
            RIDs[0].hwndTarget = Window;

            RIDs[1].usUsagePage = 0x01;
            RIDs[1].usUsage = 0x02; // NOTE(georgy): mouse input
            RIDs[1].dwFlags = 0;
            RIDs[1].hwndTarget = Window;

            if(!RegisterRawInputDevices(RIDs, ArrayCount(RIDs), sizeof(RIDs[0])))
            {
                GlobalRunning = false;
                MessageBox(Window, "Can't register input devices.", 0, MB_OK);
            }
 
            game_input GameInput = {};
            POINT MouseP;
            GetCursorPos(&MouseP);
            ScreenToClient(Window, &MouseP);
            GameInput.MouseX = MouseP.x;
            GameInput.MouseY = MouseP.y;

			while(GlobalRunning)
			{
                for(u32 ButtonIndex = 0;
                    ButtonIndex < ArrayCount(GameInput.Buttons);
                    ButtonIndex++)
                {
                    GameInput.Buttons[ButtonIndex].HalfTransitionCount = 0;
                }

				MSG Message;
				while(PeekMessage(&Message, 0, 0, 0, PM_REMOVE))
				{
                    switch(Message.message)
                    {
                        case WM_QUIT:
                        {
                            GlobalRunning = false;
                        } break;

                        case WM_INPUT:
                        {
                            RAWINPUT RawInput;
                            UINT RawInputSize = sizeof(RawInput);
                            GetRawInputData((HRAWINPUT)Message.lParam, RID_INPUT, &RawInput, &RawInputSize, sizeof(RAWINPUTHEADER));

                            if(RawInput.header.dwType == RIM_TYPEKEYBOARD)
                            {
                                UINT KeyCode = RawInput.data.keyboard.VKey;

                                b32 IsDown = (RawInput.data.keyboard.Flags == 0);
                                if(KeyCode == 'W')
                                {
                                    WinProcessKey(&GameInput.MoveForward, IsDown);
                                }
                                if(KeyCode == 'S')
                                {
                                    WinProcessKey(&GameInput.MoveBack, IsDown);
                                }
                                if(KeyCode == 'D')
                                {
                                    WinProcessKey(&GameInput.MoveRight, IsDown);
                                }
                                if(KeyCode == 'A')
                                {
                                    WinProcessKey(&GameInput.MoveLeft, IsDown);
                                }

                                if(KeyCode == VK_F4)
                                {
                                    WinProcessKey(&GameInput.F4, IsDown);
                                }
                                if(KeyCode == VK_MENU)
                                {
                                    WinProcessKey(&GameInput.Alt, IsDown);
                                }
                                if(GameInput.F4.EndedDown && GameInput.Alt.EndedDown)
                                {
                                    GlobalRunning = false;
                                }
                            }
							else if(RawInput.header.dwType == RIM_TYPEMOUSE)
                            {
                                if(RawInput.data.mouse.usFlags & MOUSE_MOVE_ABSOLUTE)
                                {
                                    GameInput.MouseXDisplacement = RawInput.data.mouse.lLastX - GameInput.MouseX;
                                    GameInput.MouseYDisplacement = RawInput.data.mouse.lLastY - GameInput.MouseY;

                                    GameInput.MouseX = RawInput.data.mouse.lLastX;
                                    GameInput.MouseY = RawInput.data.mouse.lLastY;
                                }
                                else if(RawInput.data.mouse.usFlags == MOUSE_MOVE_RELATIVE)
                                {
                                    GameInput.MouseXDisplacement = RawInput.data.mouse.lLastX;
                                    GameInput.MouseYDisplacement = RawInput.data.mouse.lLastY;

                                    GameInput.MouseX += GameInput.MouseXDisplacement;
                                    GameInput.MouseY += GameInput.MouseYDisplacement;
                                }

                                USHORT ButtonFlags = RawInput.data.mouse.usButtonFlags;
                                if((ButtonFlags & RI_MOUSE_RIGHT_BUTTON_DOWN) ||
                                   (ButtonFlags & RI_MOUSE_RIGHT_BUTTON_UP))
                                {
                                    WinProcessKey(&GameInput.MouseRight, ButtonFlags & RI_MOUSE_RIGHT_BUTTON_DOWN);
                                }
                                if((ButtonFlags & RI_MOUSE_LEFT_BUTTON_DOWN) ||
                                   (ButtonFlags & RI_MOUSE_LEFT_BUTTON_UP))
                                {
                                    WinProcessKey(&GameInput.MouseLeft, ButtonFlags & RI_MOUSE_LEFT_BUTTON_DOWN);
                                }
                            }
                        } break;

                        default:
                        {
                            TranslateMessage(&Message);
                            DispatchMessage(&Message);
                        }
                    }
				}

                RECT Rect;
                GetClientRect(Window, &Rect);
                u32 WindowWidth = Rect.right - Rect.left;
                u32 WindowHeight = Rect.bottom - Rect.top;

                glViewport(0, 0, WindowWidth, WindowHeight);
                glClearColor(1.0f, 0.0f, 0.0f, 1.0f);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                static shader Shader("shaders/Shader.glsl");

                static b32 TriDataInitialized = false;
                static GLuint TriVAO = 0, TriVBO = 0;
                static GLuint CubeVAO = 0, CubeVBO = 0;
                if(!TriDataInitialized)
                {
                    r32 TriVertices[] = 
                    {
                        -0.5f, -0.5f, 0.0f, 
                        0.5f, -0.5f, 0.0f,
                        0.0f, 0.5f, 0.0f
                    };

                    glGenVertexArrays(1, &TriVAO);
                    glGenBuffers(1, &TriVBO);
                    glBindVertexArray(TriVAO);
                    glBindBuffer(GL_ARRAY_BUFFER, TriVBO);
                    glBufferData(GL_ARRAY_BUFFER, sizeof(TriVertices), TriVertices, GL_STATIC_DRAW);
                    glEnableVertexAttribArray(0);
                    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(r32), (void *)0);
                    glBindVertexArray(0);
                    glBindBuffer(GL_ARRAY_BUFFER, 0);
                    
                    r32 CubeVertices[] = {
                        // Back face
                        -0.5f, -0.5f, -0.5f,
                        0.5f,  0.5f, -0.5f,
                        0.5f, -0.5f, -0.5f,
                        0.5f,  0.5f, -0.5f,
                        -0.5f, -0.5f, -0.5f,
                        -0.5f,  0.5f, -0.5f,
                        // Front face
                        -0.5f, -0.5f,  0.5f,
                        0.5f, -0.5f,  0.5f,
                        0.5f,  0.5f,  0.5f,
                        0.5f,  0.5f,  0.5f,
                        -0.5f,  0.5f,  0.5f,
                        -0.5f, -0.5f,  0.5f,
                        // Left face
                        -0.5f,  0.5f,  0.5f,
                        -0.5f,  0.5f, -0.5f,
                        -0.5f, -0.5f, -0.5f,
                        -0.5f, -0.5f, -0.5f,
                        -0.5f, -0.5f,  0.5f,
                        -0.5f,  0.5f,  0.5f,
                        // Right face
                        0.5f,  0.5f,  0.5f,
                        0.5f, -0.5f, -0.5f,
                        0.5f,  0.5f, -0.5f,
                        0.5f, -0.5f, -0.5f,
                        0.5f,  0.5f,  0.5f,
                        0.5f, -0.5f,  0.5f,
                        // Bottom face
                        -0.5f, -0.5f, -0.5f,
                        0.5f, -0.5f, -0.5f,
                        0.5f, -0.5f,  0.5f,
                        0.5f, -0.5f,  0.5f,
                        -0.5f, -0.5f,  0.5f,
                        -0.5f, -0.5f, -0.5f,
                        // Top face
                        -0.5f,  0.5f, -0.5f,
                        0.5f,  0.5f , 0.5f,
                        0.5f,  0.5f, -0.5f,
                        0.5f,  0.5f,  0.5f,
                        -0.5f,  0.5f, -0.5f,
                        -0.5f,  0.5f,  0.5f,
                    };
                    glGenVertexArrays(1, &CubeVAO);
                    glGenBuffers(1, &CubeVBO);
                    glBindVertexArray(CubeVAO);
                    glBindBuffer(GL_ARRAY_BUFFER, CubeVBO);
                    glBufferData(GL_ARRAY_BUFFER, sizeof(CubeVertices), CubeVertices, GL_STATIC_DRAW);
                    glEnableVertexAttribArray(0);
                    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(r32), (void *)0);
                    glBindVertexArray(0);
                    glBindBuffer(GL_ARRAY_BUFFER, 0);

                    TriDataInitialized = true;
                }

                static r32 CameraPitch = 0.0f;
                static r32 CameraHead = 0.0f;
                r32 CameraRotationSensetivity = 0.1f;
                CameraPitch -= GameInput.MouseYDisplacement*CameraRotationSensetivity;
                CameraHead -= GameInput.MouseXDisplacement*CameraRotationSensetivity;

                CameraPitch = CameraPitch > 89.0f ? 89.0f : CameraPitch;
                CameraPitch = CameraPitch < -89.0f ? -89.0f : CameraPitch;

                r32 PitchRadians = Radians(CameraPitch);
                r32 HeadRadians = Radians(CameraHead);

                r32 CameraDistanceFromHero = 5.0f;
                r32 FloorDistanceFromHero = CameraDistanceFromHero * Cos(-PitchRadians);
                
                r32 XOffsetFromHero = FloorDistanceFromHero * Sin(HeadRadians);
                r32 YOffsetFromHero = CameraDistanceFromHero * Sin(-PitchRadians);
                r32 ZOffsetFromHero = FloorDistanceFromHero * Cos(HeadRadians);
                vec3 CameraOffsetFromHero = vec3(XOffsetFromHero, YOffsetFromHero, ZOffsetFromHero);
    
                mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
                mat4 View = ViewRotationMatrixFromDirection(-CameraOffsetFromHero) * Translation(-CameraOffsetFromHero);
                mat4 Model = Identity();

                Shader.Use();
                Shader.SetMat4("Projection", Projection);
                Shader.SetMat4("View", View);
                Shader.SetMat4("Model", Model);
                glBindVertexArray(CubeVAO);
                glDrawArrays(GL_TRIANGLES, 0, 36);
                glBindVertexArray(0);
                Model = Translation(vec3(0.0f, 0.0f, -3.0f));
                Shader.SetMat4("Model", Model);
                glBindVertexArray(TriVAO);
                glDrawArrays(GL_TRIANGLES, 0, 3);
                glBindVertexArray(0);

                HDC WindowDC = GetDC(Window);
                SwapBuffers(WindowDC);
                ReleaseDC(Window, WindowDC);
			}
        }
    }

    return(0);
}