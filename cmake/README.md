OpenMS Build System Helpers
======

This folder contains all scripts used by the OpenMS build system. Please add
only build system relevant scripts. All external scripts should go into
separate directories (e.g., `cmake/modules` for `Find*` scripts).

## MSVC `import std` experiment

`OPENMS_IMPORT_STD` is OFF by default. Enabling it requires MSVC, CMake 3.30+
and a Ninja generator. It substitutes `import std;` for selected C++ standard
headers in `StringUtils.cpp`, `DateTime.cpp` and `File.cpp`. C headers needed
for macros or global C names remain included. Installed headers and pyOpenMS
are not converted. Unity builds are not supported by the experiment.

Use an x64 Visual Studio developer shell. The `MSVC import std experiment`
workflow pins CMake 3.31.10 and runs both modes on one `windows-2025` runner.
It runs for pull requests changing its workflow or `cmake/import_std.cmake`,
and can also be triggered manually. Its result artifact contains
configure/build/test logs, the original `.ninja_log` files and timings.
A failed experiment does not block the ordinary CI jobs.

For a local comparison, configure two **new** build directories with the same
preset, compiler and options, changing only `-DOPENMS_IMPORT_STD=OFF`/`ON`.
Clear `CMAKE_C_COMPILER_LAUNCHER` and `CMAKE_CXX_COMPILER_LAUNCHER` in both
configurations and disable ccache. Share a `VCPKG_INSTALLED_DIR` using the same
triplet so dependency installation can be completed before timing either build.
Then time `cmake --build <build> --target OpenMS --parallel <jobs>` with the same
job count in both directories. Do not include configuration or dependency
installation in the wall-clock measurement; include the standard library BMI
build in the ON measurement.

For each converted source, subtract the start column from the end column in
its object-file row in `.ninja_log` (milliseconds). Record the `std.ixx.obj`
and `std.compat.ixx.obj` rows separately as the one-time module compilation
cost. These per-command durations overlap in a parallel build; their sum is
not the total build wall time. Save the log before building tests, then build
`StringUtils_test`, `DateTime_test`, and `File_test` and run:

```sh
ctest --test-dir <build> -R "^(StringUtils|DateTime|File)_test$" --output-on-failure
```

Report the commit, CMake/compiler versions, configuration, job count, clean
build times, per-source and BMI times, and any compiler/linker errors in
issue #10876. Failed builds are diagnostic results, not valid speedup
measurements. Keep all options fixed between modes and repeat successful
pairs before drawing conclusions from small timing differences. This
experiment makes no promise of a speedup and does not enable other compilers.
