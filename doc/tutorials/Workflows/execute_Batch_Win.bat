@echo off
setlocal enabledelayedexpansion
 
if "%~1"=="" (
    echo Usage: %0 KNIME_HOME KNIME_WORKSPACE WORKFLOW_DIRECTORY
	echo KNIME_HOME is the directory containing knime.exe
	echo Test Data should be in KNIME_WORKSPACE\Example_Data
	echo A KNIME_WORKSPACE\TMP folder should exist for temp files
    exit /b 1
)
 
set KNIME_HOME=%~1
set WORKSPACE=%~2
set WORKFLOWS=%~3
 
if not exist "%KNIME_HOME%\knime.exe" (
    echo Error: KNIME_HOME directory "%KNIME_HOME%" does not contain knime.exe
    exit /b 1
)
 
if not exist "%WORKSPACE%" (
    echo Error: directory does not exist
    exit /b 1
)

if not exist "%WORKFLOWS%" (
    echo Error: workflow directory does not exist
    exit /b 1
)

echo KNIME_HOME = %KNIME_HOME%
echo WORKSPACE = %WORKSPACE%
echo WORKFLOWS = %WORKFLOWS%
 
 
for %%f in (%WORKFLOWS%\*.knwf) do (
    set wf=%%~nxf
	echo !wf!
	::echo %WORKFLOWS%\!wf!
    "%KNIME_HOME%\knime.exe" -consoleLog -nosplash -data "%WORKSPACE%" -application org.knime.product.KNIME_BATCH_APPLICATION --launcher.suppressErrors -nosave -workflowDir=%WORKFLOWS% -workflowFile="%WORKFLOWS%\!wf!" -reset -vmargs -Djava.io.tmpdir="%WORKSPACE%\TMP" > "%WORKSPACE%\TMP\log_!wf!.txt" 2>&1
)
 
echo Done.

::add -noexit to keep each workflow terminal alive