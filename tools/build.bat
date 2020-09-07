@echo off
REM
REM  no-fuzz Windows build script (similar to 'make' on Linux)
REM
REM  Call build.bat without arguments to print usage
REM
REM Author: Chris Bielow
REM 

IF "%~1"=="" (
  ECHO.
  ECHO.  This build script will use ALL your CPU cores ^(but on low priority^).
  ECHO.
  ECHO   Usage: build ^<target^(- for all^)^> [[^<[r]elease^|[d]ebug^|[rd]RelWithDebug^|[rm]MinSizeRel^>] ^<Sln:[a]ll^|[c]lass-test^|[t]opp^|[u]til^|[g]ui^|[d]oc^>]
  ECHO.
  ECHO  e.g.
  ECHO          // build all targets from all projects ^(TOPP, UTILS, tests, GUI^) in release mode
  ECHO          build -
  ECHO.
  ECHO          // build all targets ^(-^) from Class tests in Release mode
  ECHO          build - r c
  ECHO.
  ECHO          // just build FeatureFinderCentroided in Debug mode ^(no need to specify where to find it^)
  ECHO          build FeatureFinderCentroided d
  ECHO.
  ECHO     Targets: any Executable/library name, 'clean', '-' ^(=ALL_BUILD^)
  goto end
)

set TARGET=%~1
IF "%~1"=="-" set TARGET=ALL_BUILD

set CFG=Release
IF "%~2"=="r" set CFG=Release
IF "%~2"=="d" set CFG=Debug
IF "%~2"=="rd" set CFG=RelWithDebInfo
IF "%~2"=="rm" set CFG=MinSizeRel


set SLN=OpenMS_host.sln
IF "%~3"=="a" set SLN=OpenMS_host.sln
IF "%~3"=="c" set SLN=src\tests\class_tests\OpenMS_class_tests.sln
IF "%~3"=="t" set SLN=src\topp\openms_topp.sln
IF "%~3"=="u" set SLN=src\utils\openms_utils.sln
IF "%~3"=="g" set SLN=src\openms_gui\openms_gui.sln
IF "%~3"=="d" set SLN=doc\OpenMS_doc.sln

echo.
echo Params:
echo Target: %1
echo    Cfg: %CFG%
echo    Sln: %SLN%
echo.

if not exist %SLN% (
  ECHO.
  ECHO The .sln file '%SLN%' was not found. This script should be invoked from the root of the build tree after configuring with cmake ^(make sure you use one of the Visual Studio Generators and *not* the nmake Generator^). Change CWD and try again!
  goto end
)

REM
REM MSBuild is usually found in C:\Windows\Microsoft.NET\Framework\v4.0.30319\
REM
where /q MSBuild.exe
if not %ERRORLEVEL%==0 (
  ECHO.
  ECHO Visual Studio's 'MSBuild.exe' was not found. Please modify this build.bat to point to the correct location or make it available in %%PATH%%.
  goto end
)

set t_start=%time%
REM Invoke MSBuild.exe
REM do not use "START /WAIT /B /LOW ..." since: 
REM   - it does not allow to cancel the job on the console (it will keep running in the background)
REM   - it might trick MSBuild.exe into assuming only single-core CPU, even if /maxcpucount is specified
MSBuild.exe %SLN% /maxcpucount /target:%TARGET% /p:Configuration=%CFG%
set t_end=%time%

ECHO Time start: %t_start%
ECHO        end: %t_end%

:end
