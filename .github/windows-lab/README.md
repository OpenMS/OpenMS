# Interactive Windows test lab

Launch a fresh GitHub-hosted Windows runner to test OpenMS packages, DLL loading, and library bundling through SSH.

## Start

1. Open **Actions → Windows test lab → Run workflow**.
2. Select the branch or tag containing the workflow. The same ref is checked out for testing.
3. Keep `windows-2025` to match the current OpenMS build/release runner, or select another available image.
4. Choose a setup:
   - `runner-only` starts with the GitHub image and Visual Studio x64 tools. Use it for testing existing installers or packages.
   - `ci-dependencies` also runs the existing OpenMS dependency action (Qt, vcpkg, third-party tools, and packaging utilities). It prepares dependencies but does not configure or build OpenMS; vcpkg libraries are resolved when you configure a CMake preset. No archive/deployment SSH secrets are passed to this action.
5. Choose the interactive time limit (30 minutes by default) and leave `interactive` checked.
6. Open the **Open interactive terminal** step and copy the SSH command from its log.
7. Connect from a computer with the private key corresponding to a public SSH key on your GitHub profile. If needed, add `-i /path/to/private-key` to the command. Never copy the private key onto the runner.

The initial terminal is MSYS2 Bash **on Windows**. OpenMS requires the native Visual Studio toolchain. Enter PowerShell, which inherits the compiler and SDK environment:

```bash
pwsh -NoLogo
```

```powershell
cd $env:LAB_ROOT
Get-Command cl.exe, dumpbin.exe, cmake.exe
```

Only the person who launches the workflow can authenticate, using their GitHub SSH keys. The connection goes through the tmate relay. This provides an interactive terminal, not RDP or a persistent desktop.

## Test packages and DLLs

Download the package or installer you want to test using its URL, or install a package from the relevant package manager. For example, test a published pyOpenMS wheel in an isolated environment:

```powershell
$venv = Join-Path $env:RUNNER_TEMP 'pyopenms-test'
python -m venv $venv
& "$venv/Scripts/python.exe" -m pip install pyopenms
cd $env:RUNNER_TEMP
& "$venv/Scripts/python.exe" -c "import pyopenms; print(pyopenms.__version__)"
& "$venv/Scripts/python.exe" -m pip check
```

For a DLL or executable:

```powershell
dumpbin /DEPENDENTS C:/path/to/application.exe
dumpbin /DEPENDENTS C:/path/to/OpenMS.dll
& C:/path/to/application.exe --help
```

GitHub's image already contains many runtimes and developer tools. Check the installed bundle's startup from outside the build tree and with a minimal runtime PATH as well; a successful run in this development environment alone does not prove that all DLLs are bundled for a clean end-user machine.

## Build debugging

Select `ci-dependencies` when starting the workflow. From the OpenMS checkout, run the project's presets as needed:

```powershell
cmake --preset windows-x64-release
cmake --build --preset windows-x64-release --parallel 2
ctest --preset windows-x64-release
```

Use `windows-x64-ci` for the preset used by the build/release workflow. Initial configuration can spend substantial time building vcpkg dependencies; a full OpenMS build is not performed automatically. Keep the compiler and Release/Debug variants consistent. The runner's image version and initial tool versions are recorded in the result artifact.

## Save results and stop

Copy logs, ZIPs, wheels, or binaries into `$env:LAB_RESULTS`. For example:

```powershell
Copy-Item C:/path/to/package.zip $env:LAB_RESULTS
```

Finish from PowerShell:

```powershell
New-Item -ItemType File -Force "$env:LAB_ROOT/continue"
```

Or from Bash:

```bash
touch "$GITHUB_WORKSPACE/continue"
```

Wait for **Save test results** to finish, then download the `windows-lab-...` artifact from the run page. Artifacts are retained for seven days and follow this repository's visibility; include only files suitable for this public repository.

Closing an SSH client can leave the session running. Use the finish command or **Cancel workflow**. The interactive step has the chosen timeout, and the whole job has a three-hour ceiling including setup. Hitting a timeout can mark the run as failed; cancellation or a hard job timeout can prevent artifact upload. Every new run starts with a fresh VM.

The workflow runs only on manual dispatch, with one active session per launching user. Waiting for SSH consumes runner time. GitHub's existing Actions allowances and billing settings apply.

## Validate the environment

Uncheck `interactive` to run a short check that compiles and executes a small native Windows program and saves the tool versions and logs. This checks the runner toolchain without building OpenMS.

## Before the first merge

GitHub requires a manually dispatched workflow to exist on the repository's default branch before it can be launched. After this workflow is merged into `develop`, it appears under Actions and can be run on refs that contain it.
