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
It runs for pull requests changing its workflow, import configuration, converted
sources, discovery tests or timing analysis, and can also be triggered manually.
Superseded PR runs are cancelled. Dependencies use the same rsync-backed archive
as Windows CI, including its public HTTP restore for fork PRs; they do not consume
the repository's GitHub Actions cache quota. Compiler/build/test failures are
reported as diagnostic warnings and retained in the artifact, without preventing
the dependency action's success-only cache upload on runs with archive credentials.
The artifact contains configure/build/test logs, original `.ninja_log` files,
per-source timings and a whole-build object-duration comparison. Ordinary CI
remains responsible for required build/test checks.

For a local comparison, configure two **new** build directories with the same
preset, compiler and options, changing only `-DOPENMS_IMPORT_STD=OFF`/`ON`.
Clear `CMAKE_C_COMPILER_LAUNCHER` and `CMAKE_CXX_COMPILER_LAUNCHER` in both
configurations and disable ccache. Share a `VCPKG_INSTALLED_DIR` using the same
triplet so dependency installation can be completed before timing either build.
Do not switch an existing OFF build directory to ON: older CMake releases cache
compiler detection without the experimental gate. Use a fresh build directory
if the diagnostic reports that detection ran without import support enabled.
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
ctest --test-dir <build> -R "^(ImportStdConfiguration|(StringUtils|DateTime|File)_test)$" --output-on-failure
```

Report the commit, CMake/compiler versions, configuration, job count, clean
build times, per-source and BMI times, and any compiler/linker errors in
issue #10876. Failed builds are diagnostic results, not valid speedup
measurements. Keep all options fixed between modes and repeat successful
pairs before drawing conclusions from small timing differences. This
experiment makes no promise of a speedup and does not enable other compilers.

### Interpret the measurements

Compare all object compilations, not only the converted sources:

```sh
python tools/ci/import_std_timings.py results/ninja-OFF.log results/ninja-ON.log
```

The report separates matching objects, the three converted sources, other
matching objects and added module objects. It excludes dependency scans,
module-map generation and linking; duplicate output entries use the latest
compilation. These elapsed command-time sums overlap and are not wall-clock or
CPU-time measurements. Keep the unmodified clean-build logs for reproducibility.

In the original local clean pair, the 864 matching objects summed to 2354.106 s
OFF and 2719.793 s ON. The three converted sources accounted for just 0.116 s
of the 365.687 s increase; the other 861 objects accounted for 365.571 s.
The two new module objects added 3.513 s. Thus the 14.4% wall-clock difference
does not isolate an import-related regression: machine load, scheduling or
another build-wide difference needs investigation. These are the original
pre-rebase measurements described in PR #10881, not a new clean benchmark.
For another clean experiment, run interleaved OFF, ON, OFF, ON builds with the
same options and compare the spread as well as the totals. Keep ordinary
incremental development builds separate so collecting clean timings does not
discard their reusable objects or dependency installation.

### Choose the next sources

The original three sources remain a baseline for include/import mixing.
Their OpenMS include chains still name all four replaced C++ headers in
`StringUtils.cpp`, and four of six in `File.cpp`. In `DateTime.cpp`, they do
not name `<chrono>` or `<tuple>`, making it the clearest candidate to benefit.
This source-level check is only a lower bound: standard-library headers may
include one another too.

For a subsequent measured conversion, `DATASTRUCTURES/Date.cpp` and
`CONCEPT/UniqueIdGenerator.cpp` are stronger candidates: their OpenMS include
chains do not directly name `<chrono>`/`<tuple>` and `<chrono>`/`<random>`,
respectively. Verify the compiler's actual include trace before converting
them, retain the required C headers, and run `Date_test` and
`UniqueIdGenerator_test`. Do not expand to public headers or other compilers.
