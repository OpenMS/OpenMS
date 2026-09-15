# $Maintainer: Timo Sachsenberg $
$ErrorActionPreference = 'Stop'

$labRoot = $env:GITHUB_WORKSPACE
if (-not $labRoot -or -not $env:RUNNER_TEMP) { throw 'Run this helper through the Windows test lab workflow.' }
$resultsPath = Join-Path $env:RUNNER_TEMP 'windows-lab-results'
$scratchPath = Join-Path $env:RUNNER_TEMP 'windows-lab-scratch'
New-Item -ItemType Directory -Force -Path $resultsPath, $scratchPath | Out-Null

# Keep these available even when a later toolchain check fails.
if ($env:GITHUB_ENV) {
    "LAB_ROOT=$labRoot" | Add-Content -LiteralPath $env:GITHUB_ENV
    "LAB_RESULTS=$resultsPath" | Add-Content -LiteralPath $env:GITHUB_ENV
    "CMAKE_PREFIX_PATH=$env:CMAKE_PREFIX_PATH" | Add-Content -LiteralPath $env:GITHUB_ENV
    "SEARCH_ENGINES_DIRECTORY=$labRoot/_thirdparty" | Add-Content -LiteralPath $env:GITHUB_ENV
}

$instructions = @'
## Windows test lab

1. Open the **Open interactive terminal** step and copy its SSH command.
2. Connect using the private SSH key registered on your GitHub account.
3. The first prompt is Bash on Windows. Run `pwsh -NoLogo` for PowerShell.
4. Run `cd $env:LAB_ROOT` to enter the OpenMS checkout at the selected workflow ref.
5. Visual Studio x64 tools are configured by the workflow. Use native MSVC tools to build OpenMS.
6. Save files you want to download under `$env:LAB_RESULTS`.
7. Finish with `New-Item -ItemType File -Force "$env:LAB_ROOT/continue"`, then download the run artifact.

The selected time limit stops unattended sessions. Cancel the workflow to stop immediately; cancellation may prevent artifact upload. Each new run starts with a fresh Windows VM. If setup failed, the SSH session lets you inspect that failure.
'@
Write-Output $instructions
if ($env:GITHUB_STEP_SUMMARY) { $instructions | Add-Content -LiteralPath $env:GITHUB_STEP_SUMMARY }

$inventory = [ordered]@{
    image = $env:ImageOS
    imageVersion = $env:ImageVersion
    windows = [System.Runtime.InteropServices.RuntimeInformation]::OSDescription
    architecture = [System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture.ToString()
    powershell = $PSVersionTable.PSVersion.ToString()
    setup = $env:LAB_SETUP
    tools = [ordered]@{}
}

foreach ($tool in @('git', 'cmake', 'ninja', 'python', 'node', 'dotnet')) {
    $command = Get-Command $tool -ErrorAction SilentlyContinue
    if ($command) {
        $inventory.tools[$tool] = (& $tool --version 2>&1 | Out-String).Trim()
    } else {
        $inventory.tools[$tool] = 'Not on PATH; additional installed tools may be available through a developer shell.'
    }
}

$inventory | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $resultsPath 'environment.json')
$inventory | ConvertTo-Json -Depth 5 | Write-Output

# Compile and run a small native executable to check the Windows C++ toolchain.
if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue)) { throw 'The Visual Studio x64 compiler is missing from PATH.' }
$probeSource = Join-Path $scratchPath 'windows-lab-probe.cpp'
$probeExe = Join-Path $scratchPath 'windows-lab-probe.exe'
@'
#include <iostream>
#include <windows.h>
int main() {
    std::cout << "Windows native build OK; process " << GetCurrentProcessId() << '\n';
    return 0;
}
'@ | Set-Content -LiteralPath $probeSource
Push-Location $scratchPath
try {
    $compileOutput = & cl.exe /nologo /EHsc $probeSource "/Fe:$probeExe" 2>&1
    $compileExitCode = $LASTEXITCODE
    $compileOutput | Tee-Object -FilePath (Join-Path $resultsPath 'native-build.txt')
    if ($compileExitCode -ne 0) { throw "Native compilation failed with exit code $compileExitCode" }
    & $probeExe | Tee-Object -FilePath (Join-Path $resultsPath 'native-run.txt')
    if ($LASTEXITCODE -ne 0) { throw 'Native executable smoke check failed.' }
} finally {
    Pop-Location
}
