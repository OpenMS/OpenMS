# vcpkg Integration Verification

## Verified Findings

### Toolchain path handling on Windows

Using a relative toolchain path such as:

```powershell
-DCMAKE_TOOLCHAIN_FILE=..\vcpkg\scripts\buildsystems\vcpkg.cmake
```

can lead to a silent misconfiguration on Windows generators. In that case CMake
continues, but the vcpkg toolchain is not actually active.

Use an absolute path instead:

```powershell
-DCMAKE_TOOLCHAIN_FILE="C:\absolute\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake"
```

Indicators that the toolchain was applied correctly:

* CMake output contains `-- Running vcpkg install`
* `CMakeCache.txt` records the expected `CMAKE_TOOLCHAIN_FILE`
* Dependencies resolve from `vcpkg_installed/<triplet>`

### Windows build environment requirement

On Windows the build fails if MSVC command-line tools such as `dumpbin.exe` are
not available in `PATH`.

Use the Visual Studio Developer Command Prompt, for example:

```text
x64 Native Tools Command Prompt for VS 2022
```

### Root cause of unnecessary Qt installs

The previous `vcpkg.json` listed `qtbase` and `qtsvg` as unconditional
dependencies. vcpkg therefore resolved Qt even for CLI-only configurations where
OpenMS was configured with `-DWITH_GUI=OFF`.

That meant the dependency graph was decided before CMake reached the existing
`WITH_GUI` checks in `cmake/cmake_findExternalLibs.cmake`.

## Fix Applied

The manifest and top-level CMake configuration now split GUI dependencies from
core dependencies:

* `vcpkg.json` declares `gui` as the manifest default feature for normal builds
* `vcpkg.json` keeps only core packages in the default dependency set
* Qt packages moved into a new manifest feature: `gui`
* `CMakeLists.txt` now defines `WITH_GUI` before the first `project()` call
* When `WITH_GUI=ON`, OpenMS automatically requests the vcpkg manifest feature
  `gui`
* When `WITH_GUI=OFF`, OpenMS requests no GUI manifest feature and sets
  `VCPKG_MANIFEST_NO_DEFAULT_FEATURES=ON`, so Qt is not part of the dependency
  solve

This preserves the default GUI-enabled behavior while allowing a true CLI-only
manifest install.

If vcpkg is invoked manually instead of through CMake, GUI builds must request
the `gui` manifest feature explicitly, and CLI-only builds must disable the
manifest default features.

## Validation Notes

Expected behavior after this change:

* Default configure path keeps GUI builds working and still pulls Qt through the
  `gui` manifest default feature
* `-DWITH_GUI=OFF` sets `VCPKG_MANIFEST_NO_DEFAULT_FEATURES=ON`
* `-DWITH_GUI=OFF` no longer requests `qtbase` or `qtsvg`
* CLI-only environments avoid the large Qt dependency build, reducing install
  time and resource usage

When switching between GUI and CLI-only builds, prefer a fresh build directory
so the CMake cache and manifest-installed dependency set stay in sync.

## Recommended Checks

For GUI-enabled validation:

```powershell
cmake -S . -B build-gui `
  -DCMAKE_TOOLCHAIN_FILE="C:\absolute\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake" `
  -DWITH_GUI=ON
```

For CLI-only validation:

```powershell
cmake -S . -B build-cli `
  -DCMAKE_TOOLCHAIN_FILE="C:\absolute\path\to\vcpkg\scripts\buildsystems\vcpkg.cmake" `
  -DWITH_GUI=OFF
```

The CLI-only configure should not install or reference Qt packages.



