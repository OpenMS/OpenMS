# pyOpenMS Build System Modernization

## Overview

This document tracks the modernization of the pyOpenMS build system from CMake-based Python dependency management to a modern Python packaging approach using `py-build-cmake` and `uv`.

## Goals

1. Remove Python dependency checks and scripts from CMake
2. Use `py-build-cmake` as build backend with `uv`
3. Add CMake config settings to `pyproject.toml`
4. Simplify CI to use `uv build --wheel`
5. Move shared lib fixup to `delvewheel`/`delocate`/`auditwheel`
6. Use `cibuildwheel` for nightly/release wheels
7. Consolidate CI jobs and add pyOpenMS to `openms-ci-full`

---

## Implementation Status

### Completed

- [x] **pyproject.toml created** (`src/pyOpenMS/pyproject.toml`)
  - py-build-cmake as build backend
  - Build dependencies: autowrap, Cython, numpy
  - Runtime dependencies: numpy, pandas, matplotlib
  - cibuildwheel configuration for all platforms
  - pytest configuration

- [x] **CMakeLists.txt modernized** (`src/pyOpenMS/CMakeLists.txt`)
  - Removed `requirements_check.py` execution
  - Simplified env.py generation (inline file write)
  - Added `install()` commands for py-build-cmake
  - Changed `NO_DEPENDENCIES` default to `ON`
  - Kept backward compatibility with legacy targets

- [x] **New CI workflow created** (`.github/workflows/pyopenms-wheels-cibuildwheel.yml`)
  - Uses cibuildwheel for multi-platform wheel building
  - Separate OpenMS build jobs per platform
  - Uses auditwheel/delocate/delvewheel for dependency bundling
  - Test job validates wheels on all platforms
  - Publish job for PyPI upload

- [x] **Quick test added to openms-ci-full** (`.github/workflows/openms_ci_matrix_full.yml`)
  - New `pyopenms-quick-test` job
  - Runs on PRs to develop
  - Tests new py-build-cmake approach with `uv build --wheel`

- [x] **Local build tested successfully** (2024-12-19)
  - Fixed Eigen3 version range in `cmake/OpenMSConfig.cmake.in`
  - Fixed share directory copy when path doesn't exist
  - Built wheel with `uv build --wheel`
  - Successfully imported and tested pyopenms

- [x] **Deprecated files removed** (2024-12-19)
  - Removed `requirements_bld.txt`, `requirements_check.py`, `setup.py`, `env.py.in`
  - Removed old `pyopenms-wheels.yml` workflow

### Pending

- [ ] Full CI validation on all platforms

---

## Files Changed

### New Files

| File | Purpose |
|------|---------|
| `src/pyOpenMS/pyproject.toml` | PEP 517/518/621 build configuration |
| `.github/workflows/pyopenms-wheels-cibuildwheel.yml` | New cibuildwheel-based CI |

### Modified Files

| File | Changes |
|------|---------|
| `src/pyOpenMS/CMakeLists.txt` | Simplified, added install() targets, removed requirements_check, added existence check for share dir |
| `.github/workflows/openms_ci_matrix_full.yml` | Added pyopenms-quick-test job |
| `cmake/OpenMSConfig.cmake.in` | Fixed Eigen3 version range handling for older CMake/Eigen configs |

### Files Removed

| File | Reason |
|------|--------|
| `src/pyOpenMS/requirements_bld.txt` | Moved to pyproject.toml |
| `src/pyOpenMS/requirements_check.py` | No longer needed |
| `src/pyOpenMS/setup.py` | Replaced by pyproject.toml |
| `src/pyOpenMS/env.py.in` | Simplified to inline generation |
| `.github/workflows/pyopenms-wheels.yml` | Replaced by pyopenms-wheels-cibuildwheel.yml |

---

## Build Instructions

### Local Development Build

```bash
# 1. Build OpenMS first
cd /path/to/OpenMS
mkdir build && cd build
cmake .. -DPYOPENMS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build . --target OpenMS OpenSwathAlgo -j4

# 2. Build pyOpenMS wheel
cd /path/to/OpenMS/src/pyOpenMS
export OpenMS_ROOT=/path/to/OpenMS/build
export CMAKE_PREFIX_PATH=/path/to/contrib:/path/to/OpenMS/build
pip install uv
uv build --wheel

# 3. Install and test
pip install dist/*.whl
python -c "import pyopenms; print(pyopenms.EmpiricalFormula('C6H12O6').getMonoWeight())"
```

### Using pip directly

```bash
cd /path/to/OpenMS/src/pyOpenMS
pip wheel . --no-build-isolation \
    --config-settings=cmake.options.OpenMS_ROOT=/path/to/OpenMS/build
```

---

## CI Workflow Details

### pyopenms-wheels-cibuildwheel.yml

```
Trigger: Release tags, nightly branch, develop PRs

Jobs:
  build-openms-linux      -> Builds OpenMS on Linux x64
  build-openms-linux-arm  -> Builds OpenMS on Linux ARM64
  build-openms-macos      -> Builds OpenMS on macOS ARM64
  build-openms-windows    -> Builds OpenMS on Windows x64

  build-wheels-linux      -> cibuildwheel + auditwheel
  build-wheels-linux-arm  -> cibuildwheel + auditwheel
  build-wheels-macos      -> cibuildwheel + delocate
  build-wheels-windows    -> cibuildwheel + delvewheel

  test-wheels             -> Install and test on all platforms
  publish                 -> Upload to PyPI
```

### openms_ci_matrix_full.yml (Addition)

```
New Job: pyopenms-quick-test
  Runs on: PRs to develop
  Container: manylinux_2_34
  Steps:
    1. Build OpenMS with PYOPENMS=OFF
    2. Build wheel with uv build --wheel
    3. Test import and basic functionality
```

---

## Configuration Reference

### pyproject.toml Key Settings

```toml
[build-system]
requires = [
    "py-build-cmake>=0.3.0",
    "autowrap>=0.24.0",
    "Cython>=3.1.0",
    "numpy>=2.0",
]
build-backend = "py_build_cmake.build"

[tool.py-build-cmake.cmake.options]
PY_NUM_MODULES = "8"
PY_NUM_THREADS = "4"
NO_DEPENDENCIES = "ON"

[tool.cibuildwheel]
build = "cp310-* cp311-* cp312-* cp313-* cp314-*"
skip = "*-win32 *-musllinux_* pp*"
```

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `OpenMS_ROOT` | Path to OpenMS build/install directory |
| `CMAKE_PREFIX_PATH` | Additional paths for dependency discovery |
| `PY_NUM_MODULES` | Number of Cython modules (default: 8) |
| `PY_NUM_THREADS` | Build parallelism (default: 4) |

---

## Wheel Repair Tools

| Platform | Tool | Command |
|----------|------|---------|
| Linux | auditwheel | `auditwheel repair -w {dest_dir} {wheel}` |
| macOS | delocate | `delocate-wheel --require-archs {arch} -w {dest_dir} {wheel}` |
| Windows | delvewheel | `delvewheel repair -w {dest_dir} {wheel}` |

---

## Progress Log

### 2024-12-19: Implementation Complete

**Completed:**
1. Created `pyproject.toml` with:
   - py-build-cmake build backend configuration
   - Full project metadata (PEP 621)
   - cibuildwheel configuration for all platforms
   - Dependency groups for build/dev/test

2. Modified `CMakeLists.txt`:
   - Removed `requirements_check.py` execution
   - Simplified env.py generation to inline file write
   - Added `install()` targets for py-build-cmake
   - Changed `NO_DEPENDENCIES` default to `ON`
   - Kept legacy `pyopenms_wheel` target for backward compatibility
   - Reduced file from 640 to 522 lines

3. Created new CI workflow:
   - `pyopenms-wheels-cibuildwheel.yml` using cibuildwheel
   - Separate OpenMS build jobs for each platform
   - Platform-specific wheel repair (auditwheel/delocate/delvewheel)
   - Test and publish jobs

4. Updated `openms_ci_matrix_full.yml`:
   - Added `pyopenms-quick-test` job for PR validation
   - Tests new build system with `uv build --wheel`

### 2024-12-19 (continued): Local Build Tested Successfully

**Fixed issues during testing:**
1. Fixed Eigen3 version range in `cmake/OpenMSConfig.cmake.in`:
   - Changed `find_dependency` to `find_package` with QUIET for version range check
   - This allows fallback when Eigen3Config.cmake doesn't support version ranges

2. Fixed share directory copy in `src/pyOpenMS/CMakeLists.txt`:
   - Added `EXISTS` check before `file(COPY)` to handle missing share directory
   - Prevents build failure when OPENMS_SHARE_DIR points to non-existent path

**Test results:**
```bash
$ uv build --wheel
# Successfully built pyopenms-3.4.0-cp310-cp310-linux_x86_64.whl (9.7 MB)

$ python3 -c "import pyopenms; print(pyopenms.EmpiricalFormula('C6H12O6').getMonoWeight())"
180.0633903828
```

### Next Steps

1. Run CI to validate all platforms (create PR to trigger workflows)

---

## Notes

### Backward Compatibility

The new CMakeLists.txt maintains backward compatibility:
- Legacy `pyopenms_wheel` target still works
- Can still use old workflow until new one is validated
- env.py is still generated for create_cpp_extension.py

### Key Technical Decisions

1. **NO_DEPENDENCIES=ON by default**: Wheel repair tools handle dependency bundling
2. **Separate OpenMS build**: OpenMS must be built before pyOpenMS
3. **env.py simplified**: Only variables needed by create_cpp_extension.py
4. **install() targets added**: Required for py-build-cmake packaging

### Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| py-build-cmake limitations | Kept legacy targets as fallback |
| cibuildwheel complexity | Pre-build OpenMS in separate jobs |
| Platform-specific issues | Thorough testing on all platforms |
| Backward compatibility | Gradual migration, old workflow kept |
