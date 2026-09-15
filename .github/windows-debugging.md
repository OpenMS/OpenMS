# Debugging Windows builds with tmate

The existing **Build and Test**, **Release**, and
**pyopenms-wheels-cibuildwheel** workflows have an optional `debug_windows`
checkbox under **Actions → Run workflow**. It is off by default.

Enable it when manually launching a build to pause the Windows runner after:

- **Build and Test:** the OpenMS build and test attempt.
- **Release:** the OpenMS build and packaging attempt.
- **pyopenms-wheels-cibuildwheel:** the Windows wheel build/repair attempt.

The session opens on success or failure, so the checkout, dependencies, build
tree, and any generated packages remain available for inspection. Linux and
macOS jobs continue normally. Pull requests and automatic pushes never open
tmate sessions. Debug runs skip installer/documentation deployment and wheel
publishing; artifacts remain available through the workflows' usual uploads
when the relevant preceding steps succeed.

## Connect

Open the Windows job's **Debug Windows build with tmate** step and copy the SSH
command from its log. Connect using the private key corresponding to a public
SSH key registered on the GitHub account that launched the workflow. If needed,
add `-i /path/to/private-key` to the command. Do not copy private keys to the runner.

The initial prompt is MSYS2 Bash on the Windows runner. Start native PowerShell:

```bash
pwsh -NoLogo
```

```powershell
cd $env:GITHUB_WORKSPACE
# Native build/release workflow:
cmake --build build/windows-x64-ci --target OpenMS
# Inspect imported DLLs in the binary or wheel staging directory you are testing:
dumpbin /DEPENDENTS path/to/OpenMS.dll
```

The Visual Studio environment is inherited when setup reached that step.
The native workflows use `build/windows-x64-ci`; the wheel workflow uses
`build`, `install`, and `wheelhouse`. A failure may leave only some of these
directories. A manual repair does not erase the original failed step's status.

## Finish

Continue the workflow from PowerShell:

```powershell
New-Item -ItemType File -Force "$env:GITHUB_WORKSPACE/continue"
```

Or from Bash: `touch "$GITHUB_WORKSPACE/continue"`.

The tmate step is limited to 30 minutes, within the job's overall time limit.
Use **Cancel workflow** to stop immediately. Closing the SSH client can leave
the session running. A timeout can fail the step, and cancellation can prevent
subsequent artifact uploads. The VM is discarded when the job ends; save any
needed files through the workflow's artifact paths before continuing.

Only keys registered to the launching GitHub user are accepted. The session
uses the third-party tmate relay and consumes normal Actions runner time.
