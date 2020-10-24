@echo off

REM For Visual Studio Code
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" x64
c:;cd Users\georg\source\repos\Haha3D\Haha3D

set CommonCompilerFlags=-MT -nologo -Gm- -GR- -EHsc -fp:fast -O2 -Oi -WX -W4 -wd4324 -wd4505 -wd4456 -wd4457 -wd4063 -wd4702 -wd4201 -wd4100 -wd4189 -wd4459 -wd4127 -wd4311 -wd4302 -FC -Z7 
set CommonCompilerFlags=-DGLEW_STATIC %CommonCompilerFlags%
set CommonCompilerFlags=-D_CRT_SECURE_NO_WARNINGS %CommonCompilerFlags%
set CommonLinkerFlags=-incremental:no -opt:ref user32.lib gdi32.lib Winmm.lib opengl32.lib glew32s.lib -ignore:4099

IF NOT EXIST ..\build mkdir ..\build
pushd ..\build

del *.pdb >NUL 2> NUL

REM 64-bit build
cl %CommonCompilerFlags% "..\Haha3D\code\win_haha3d.cpp" -Fmwin_haha3d.map /link -PDB:win_haha3d.pdb %CommonLinkerFlags%

popd