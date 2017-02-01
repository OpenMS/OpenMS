@echo off
REM
REM  no-fuzz Windows build script (similar to 'make' on Linux)
REM
REM
REM Author: Chris Bielow
REM 

IF "%~1"=="" (
  ECHO.
  ECHO.  This build script will use ALL your CPU cores ^(but on low priority^).
  ECHO.
  ECHO   Usage: build ^<target^(- for all^)^> [[^<[r]elease^|[d]ebug^>] ^<Sln:[a]ll^|[c]lass-test^|[t]opp^|[u]til^|[g]ui^>]
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
  goto end
)

set TARGET=%~1
IF "%~1"=="-" set TARGET=ALL_BUILD


IF "%~2"=="" set CFG=Release
IF "%~2"=="r" set CFG=Release
IF "%~2"=="d" set CFG=Debug


IF "%~3"==""  set SLN=OpenMS_host.sln
IF "%~3"=="a" set SLN=OpenMS_host.sln
IF "%~3"=="c" set SLN=src\tests\class_tests\OpenMS_class_tests.sln
IF "%~3"=="t" set SLN=src\topp\openms_topp.sln
IF "%~3"=="u" set SLN=src\utils\openms_utils.sln
IF "%~3"=="g" set SLN=src\openms_gui\openms_gui.sln

echo.
echo Params:
echo Target: %1
echo    Cfg: %CFG%
echo    Sln: %SLN%
echo.

if not exist %SLN% (
  ECHO.
  ECHO The .sln file '%%SLN%%' was not found. This script should be invoked from the root of the build tree. Change CWD and try again!
  goto end
)

REM
REM MSBuild is usually found in C:\Windows\Microsoft.NET\Framework\v4.0.30319\
REM
where /q MSBuild.exe
if not %ERRORLEVEL%==0 (
  ECHO.
  ECHO Visual Studio's 'MSBuild.exe' was not found. Please modify this .bat file to point to the correct location or make it available in %%PATH%%.
  goto end
)
REM run with low priority, so the machine is usable, but on full steam when idle
start /WAIT /B /LOW MSBuild.exe %SLN% /maxcpucount /target:%TARGET% /p:Configuration=%CFG%


:end
