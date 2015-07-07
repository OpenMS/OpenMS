@echo off
REM
REM  no-fuzz Windows build script (similar to 'make' on Linux)
REM
REM
REM Author: Chris Bielow
REM 

IF "%~1"=="" (
  ECHO.
  ECHO   Usage: b.bat ^<target^(- for all^)^> [[^<[r]elease^|[d]ebug^>] ^<Sln:[c]lass-test^|[t]opp^|[u]til^|[g]ui^>]
  ECHO.
  ECHO  e.g.
  ECHO          b.bat -    
  ECHO          b.bat - r c   // same as above ^(build all class tests in release mode^)
  ECHO          b.bat FeatureFinderCentroided r t
  goto end
)

set TARGET=%~1
IF "%~1"=="-" set TARGET=ALL_BUILD


IF "%~2"=="" set CFG=Release
IF "%~2"=="r" set CFG=Release
IF "%~2"=="d" set CFG=Debug


IF "%~3"=="" set SLN=src\tests\class_tests\OpenMS_class_tests.sln
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
  ECHO The .sln file '%SLN%' was not found. This script should be invoked from the root of the build tree. Change CWD and try again!
  goto end
)

REM
REM MSBUild is usually found in C:\Windows\Microsoft.NET\Framework\v4.0.30319\
REM
where /q MSBuild.exe
if not %ERRORLEVEL%==0 (
  ECHO.
  ECHO Visual Studio's 'MSBuild.exe' was not found. Please modify this .bat file to point to the correct location or make it available in $PATH.
  goto end
)
MSBuild.exe %SLN% /target:%TARGET% /p:Configuration=%CFG%


:end