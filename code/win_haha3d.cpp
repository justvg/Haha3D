#include <Windows.h>

#include <string>
#include <fstream>

#include "glew\glew.h"
#include "glew\wglew.h"

#include "haha3d_platform.h"
#include "haha3d_renderer_opengl.cpp"

global_variable b32 GlobalRunning;
global_variable b32 GlobalWindowIsFocused;
global_variable LARGE_INTEGER GlobalPerformanceFrequency;

internal void
CatStrings(u32 ALength, char *SourceA, u32 BLength, char *SourceB, u32 DestSize, char *Dest)
{
    Assert(DestSize > (ALength + BLength));

    for(u32 Index = 0;
        Index < ALength;
        Index++)
    {
        *Dest++ = *SourceA++;
    }

    for(u32 Index = 0;
        Index < BLength;
        Index++)
    {
        *Dest++ = *SourceB++;
    }

    *Dest = 0;
}

struct win_game_code
{
    HMODULE GameCodeDLL;
    FILETIME DLLLastWriteTime;

    game_update_and_render *UpdateAndRender;
};

internal FILETIME 
WinGetLastWriteTime(char *Filename)
{
    FILETIME LastWriteTime = {};

    WIN32_FILE_ATTRIBUTE_DATA Data;
    if(GetFileAttributesEx(Filename, GetFileExInfoStandard, &Data))
    {
        LastWriteTime = Data.ftLastWriteTime;
    }

    return(LastWriteTime);
}

internal win_game_code
WinLoadGameCode(char *SourceDLLName, char *TempDLLName)
{
    win_game_code Result = {};

    Result.DLLLastWriteTime = WinGetLastWriteTime(SourceDLLName);

    CopyFile(SourceDLLName, TempDLLName, FALSE);

    Result.GameCodeDLL = LoadLibrary(TempDLLName);
    if(Result.GameCodeDLL)
    {
        Result.UpdateAndRender = (game_update_and_render *)GetProcAddress(Result.GameCodeDLL, "GameUpdateAndRender");
    }

    return(Result);
}

internal void
WinUnloadGameCode(win_game_code *GameCode)
{
    if(GameCode->GameCodeDLL)
    {
        FreeLibrary(GameCode->GameCodeDLL);
        GameCode->GameCodeDLL = 0;
    }

    GameCode->UpdateAndRender = 0;
}

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

            wglSwapIntervalEXT(1);
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
                GlobalWindowIsFocused = true;
            }
            else
            {
                ShowCursor(TRUE);
                GlobalWindowIsFocused = false;
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

inline void
WinProcessKey(button *Button, b32 IsDown)
{
    if(Button->EndedDown != IsDown)
    {  
        Button->HalfTransitionCount++;
    }

    Button->EndedDown = IsDown;
}

inline LARGE_INTEGER
WinGetPerformanceCounter(void)
{
    LARGE_INTEGER Result;
    QueryPerformanceCounter(&Result);
    return(Result);
}

inline r32
WinGetSecondsElapsed(LARGE_INTEGER Start, LARGE_INTEGER End)
{
    r32 Result = (End.QuadPart - Start.QuadPart) / (r32)GlobalPerformanceFrequency.QuadPart;
    return(Result);
}

int CALLBACK
WinMain(HINSTANCE Instance, HINSTANCE PrevInstance, LPSTR CommandLine, int ShowCode)
{
    QueryPerformanceCounter(&GlobalPerformanceFrequency);

    UINT DesiredSchedulerMS = 1;
    b32 SleepIsGranular = (timeBeginPeriod(DesiredSchedulerMS) == TIMERR_NOERROR);

    char ExecutableFilename[MAX_PATH];
    char *OneAfterLastExecutableFilenameSlash = 0;
    DWORD ExecutableFilenameLength = GetModuleFileName(0, ExecutableFilename, sizeof(ExecutableFilename));
    for(char *Scan = ExecutableFilename;
        *Scan;
        Scan++)
    {
        if(*Scan == '\\')
        {
            OneAfterLastExecutableFilenameSlash = Scan + 1;
        }
    }

    char GameCodeSourceDLLName[MAX_PATH];
    char GameCodeTempDLLName[MAX_PATH];
    CatStrings((u32)(OneAfterLastExecutableFilenameSlash - ExecutableFilename), ExecutableFilename,
                sizeof("haha3d.dll"), "haha3d.dll", 
                sizeof(GameCodeSourceDLLName), GameCodeSourceDLLName);
    CatStrings((u32)(OneAfterLastExecutableFilenameSlash - ExecutableFilename), ExecutableFilename,
                sizeof("haha3d_temp.dll"), "haha3d_temp.dll", 
                sizeof(GameCodeTempDLLName), GameCodeTempDLLName);

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

            platform_api PlatformAPI = {};
            PlatformAPI.InitBuffers = InitBuffers;
            PlatformAPI.InitBuffersWithEBO = InitBuffersWithEBO;
            PlatformAPI.CompileShader = CompileShader;

            game_memory GameMemory = {};

            GameMemory.PermanentStorageSize = Megabytes(256);
            GameMemory.PermanentStorage = VirtualAlloc(0, GameMemory.PermanentStorageSize, MEM_COMMIT|MEM_RESERVE, PAGE_READWRITE);

            GameMemory.PlatformAPI = PlatformAPI;

            win_game_code Game = WinLoadGameCode(GameCodeSourceDLLName, GameCodeTempDLLName);
            
            r32 TargetSecondsPerFrame = 1.0f / 60.0f;
            GameInput.dtForFrame = TargetSecondsPerFrame;
            LARGE_INTEGER LastCounter = WinGetPerformanceCounter();
			while(GlobalRunning)
			{
                if(GlobalWindowIsFocused)
                {
                    RECT ClipRect;
                    GetWindowRect(Window, &ClipRect);
                    LONG ClipRectHeight = ClipRect.bottom - ClipRect.top;
                    LONG ClipRectWidth = ClipRect.right - ClipRect.left;
                    ClipRect.top = ClipRect.bottom = ClipRect.top + ClipRectHeight/2;
                    ClipRect.left = ClipRect.right = ClipRect.left + ClipRectWidth/2;
                    ClipCursor(&ClipRect);
                }
                else
                {
                    ClipCursor(0);
                }

                GameInput.MouseXDisplacement = GameInput.MouseYDisplacement = 0;
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

                render_command_buffer RenderCommandBuffer;
                RenderCommandBuffer.CommandCount = 0;
                RenderCommandBuffer.CurrentShaderID = 0;
                RenderCommandBuffer.ShadersStackTop = 0;

                if(Game.UpdateAndRender)
                {
                    Game.UpdateAndRender(&GameMemory, &GameInput, &RenderCommandBuffer, WindowWidth, WindowHeight);
                }

                HDC WindowDC = GetDC(Window);
                RenderCommands(&RenderCommandBuffer);
                SwapBuffers(WindowDC);
                ReleaseDC(Window, WindowDC);

                FILETIME NewGameCodeDLLWriteTime = WinGetLastWriteTime(GameCodeSourceDLLName);
                b32 GameCodeNeedsToBeReloaded = CompareFileTime(&NewGameCodeDLLWriteTime, &Game.DLLLastWriteTime);
                if(GameCodeNeedsToBeReloaded)
                {
                    WinUnloadGameCode(&Game);
                    Game = WinLoadGameCode(GameCodeSourceDLLName, GameCodeTempDLLName);
                }

#if 0
                r32 SecondsElapsedForFrame = WinGetSecondsElapsed(LastCounter, WinGetPerformanceCounter());
                if(SecondsElapsedForFrame < TargetSecondsPerFrame)
                {
                    if(SleepIsGranular)
                    {
                        DWORD SleepMS = (DWORD)(1000.0f * (TargetSecondsPerFrame - SecondsElapsedForFrame));

                        if(SleepMS > 0)
                        {
                            Sleep(SleepMS);
                        }
                    }

                    while(SecondsElapsedForFrame < TargetSecondsPerFrame)
                    {
                        SecondsElapsedForFrame = WinGetSecondsElapsed(LastCounter, WinGetPerformanceCounter());
                    }
                }
#endif
                LastCounter = WinGetPerformanceCounter();
			}
        }
    }

    return(0);
}
