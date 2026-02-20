# OpenMS Library Decomposition Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Split the monolithic `libOpenMS.so` into three sub-libraries (`OpenMS_Core`, `OpenMS_IO`, `OpenMS_Algo`) for faster incremental builds and enforced dependency boundaries.

**Architecture:** Three libraries with strict downward-only layering: `Algo → IO → Core`. A compatibility `OpenMS` INTERFACE target links all three transitively so existing consumers don't break. Source files and headers stay in their current directories — only CMake wiring changes.

**Tech Stack:** CMake 3.15+, C++20, `openms_add_library()` macro from `cmake/add_library_macros.cmake`, CTest

**Design doc:** `docs/plans/2026-02-20-library-decomposition-design.md`

---

## Phase 1: Establish OpenMS_Core Library

### Task 1: Create Core includes.cmake

**Files:**
- Create: `src/openms/core_includes.cmake`

**Context:** The current `src/openms/includes.cmake` aggregates all 67 `sources.cmake` files into `OpenMS_sources` and `OpenMS_sources_h`. We need a separate aggregation file for Core modules only.

**Step 1: Create core_includes.cmake**

Core modules: CONCEPT, SYSTEM, MATH, DATASTRUCTURES, METADATA, KERNEL, CHEMISTRY, IONMOBILITY

```cmake
# core_includes.cmake — aggregates sources for OpenMS_Core library
set(OpenMS_Core_sources CACHE INTERNAL "OpenMS Core source files")

## Source files — order respects dependency hierarchy
include(source/INTERFACES_IMPL/sources.cmake)
include(source/CONCEPT/sources.cmake)
include(source/SYSTEM/sources.cmake)
include(source/MATH/sources.cmake)
include(source/MATH/STATISTICS/sources.cmake)
include(source/MATH/MISC/sources.cmake)
include(source/DATASTRUCTURES/sources.cmake)
include(source/METADATA/sources.cmake)
include(source/METADATA/ID/sources.cmake)
include(source/KERNEL/sources.cmake)
include(source/CHEMISTRY/sources.cmake)
include(source/CHEMISTRY/ISOTOPEDISTRIBUTION/sources.cmake)
include(source/CHEMISTRY/MASSDECOMPOSITION/IMS/sources.cmake)
include(source/CHEMISTRY/MASSDECOMPOSITION/sources.cmake)
include(source/IONMOBILITY/sources.cmake)

# Rename the accumulated variable from OpenMS_sources to OpenMS_Core_sources
set(OpenMS_Core_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS Core source files")
# Reset OpenMS_sources so the next includes.cmake starts fresh
set(OpenMS_sources CACHE INTERNAL "")

## Header files
set(OpenMS_Core_sources_h CACHE INTERNAL "OpenMS Core header files")

include(include/OpenMS/INTERFACES/sources.cmake)
include(include/OpenMS/CONCEPT/sources.cmake)
include(include/OpenMS/SYSTEM/sources.cmake)
include(include/OpenMS/MATH/MISC/sources.cmake)
include(include/OpenMS/MATH/sources.cmake)
include(include/OpenMS/MATH/STATISTICS/sources.cmake)
include(include/OpenMS/DATASTRUCTURES/sources.cmake)
include(include/OpenMS/DATASTRUCTURES/Utils/sources.cmake)
include(include/OpenMS/METADATA/sources.cmake)
include(include/OpenMS/METADATA/ID/sources.cmake)
include(include/OpenMS/KERNEL/sources.cmake)
include(include/OpenMS/CHEMISTRY/sources.cmake)
include(include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/sources.cmake)
include(include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/sources.cmake)
include(include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/sources.cmake)
include(include/OpenMS/IONMOBILITY/sources.cmake)

set(OpenMS_Core_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS Core header files")
set(OpenMS_sources_h CACHE INTERNAL "")
```

**Important note:** The existing `sources.cmake` files append to `OpenMS_sources` and `OpenMS_sources_h`. We need to either:
- (a) Modify each `sources.cmake` to accept a variable name parameter, or
- (b) Use the reset-and-capture pattern shown above, or
- (c) Refactor `sources.cmake` files to use a configurable list name

Option (b) is simplest for Phase 1. However, this requires that each `includes.cmake` (core, io, algo) is included sequentially and resets the accumulator variables between uses.

**Step 2: Validate the source list is complete**

Run a quick check that all expected modules are included. Count the files:
```bash
# From src/openms/
find source/KERNEL source/CONCEPT source/DATASTRUCTURES source/MATH source/SYSTEM source/METADATA source/CHEMISTRY source/IONMOBILITY source/INTERFACES_IMPL -name "*.cpp" | wc -l
```
Expected: ~240 .cpp files.

**Step 3: Commit**

```bash
git add src/openms/core_includes.cmake
git commit -m "build: add core_includes.cmake for OpenMS_Core library sources"
```

---

### Task 2: Create IO includes.cmake

**Files:**
- Create: `src/openms/io_includes.cmake`

**Context:** IO library = FORMAT module and its subdirectories. Specialized format handlers that depend on ANALYSIS will move to Algo in a later task; for now, include all of FORMAT.

**Step 1: Create io_includes.cmake**

```cmake
# io_includes.cmake — aggregates sources for OpenMS_IO library
set(OpenMS_IO_sources CACHE INTERNAL "OpenMS IO source files")

## Source files
include(source/FORMAT/DATAACCESS/sources.cmake)
include(source/FORMAT/HANDLERS/sources.cmake)
include(source/FORMAT/MSNUMPRESS/sources.cmake)
include(source/FORMAT/VALIDATORS/sources.cmake)
include(source/FORMAT/OPTIONS/sources.cmake)
include(source/FORMAT/sources.cmake)

set(OpenMS_IO_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS IO source files")
set(OpenMS_sources CACHE INTERNAL "")

## Header files
set(OpenMS_IO_sources_h CACHE INTERNAL "OpenMS IO header files")

include(include/OpenMS/FORMAT/sources.cmake)
include(include/OpenMS/FORMAT/DATAACCESS/sources.cmake)
include(include/OpenMS/FORMAT/HANDLERS/sources.cmake)
include(include/OpenMS/FORMAT/MSNUMPRESS/sources.cmake)
include(include/OpenMS/FORMAT/VALIDATORS/sources.cmake)
include(include/OpenMS/FORMAT/OPTIONS/sources.cmake)

set(OpenMS_IO_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS IO header files")
set(OpenMS_sources_h CACHE INTERNAL "")
```

**Step 2: Commit**

```bash
git add src/openms/io_includes.cmake
git commit -m "build: add io_includes.cmake for OpenMS_IO library sources"
```

---

### Task 3: Create Algo includes.cmake

**Files:**
- Create: `src/openms/algo_includes.cmake`

**Context:** Algo = everything not in Core or IO: ML, ANALYSIS, PROCESSING, FEATUREFINDER, COMPARISON, QC, APPLICATIONS, and OPENSWATH (conditional).

**Step 1: Create algo_includes.cmake**

```cmake
# algo_includes.cmake — aggregates sources for OpenMS_Algo library
set(OpenMS_Algo_sources CACHE INTERNAL "OpenMS Algo source files")

## Source files
include(source/ML/sources.cmake)
include(source/ML/NNLS/sources.cmake)
include(source/ML/SVM/sources.cmake)
include(source/ML/CLUSTERING/sources.cmake)
include(source/ML/GRIDSEARCH/sources.cmake)
include(source/ML/CROSSVALIDATION/sources.cmake)
include(source/ML/INTERPOLATION/sources.cmake)
include(source/ML/ROCCURVE/sources.cmake)
include(source/ML/RANSAC/sources.cmake)
include(source/ML/REGRESSION/sources.cmake)
include(source/ANALYSIS/QUANTITATION/sources.cmake)
include(source/ANALYSIS/SEQUENCE/sources.cmake)
include(source/ANALYSIS/MAPMATCHING/sources.cmake)
include(source/ANALYSIS/DECHARGING/sources.cmake)
include(source/ANALYSIS/ID/sources.cmake)
include(source/ANALYSIS/MRM/sources.cmake)
include(source/ANALYSIS/NUXL/sources.cmake)
include(source/ANALYSIS/TARGETED/sources.cmake)
include(source/ANALYSIS/TOPDOWN/sources.cmake)
include(source/ANALYSIS/XLMS/sources.cmake)
include(source/PROCESSING/BASELINE/sources.cmake)
include(source/PROCESSING/FILTERING/sources.cmake)
include(source/PROCESSING/RESAMPLING/sources.cmake)
include(source/PROCESSING/MISC/sources.cmake)
include(source/PROCESSING/DEISOTOPING/sources.cmake)
include(source/PROCESSING/CALIBRATION/sources.cmake)
include(source/PROCESSING/FEATURE/sources.cmake)
include(source/PROCESSING/SMOOTHING/sources.cmake)
include(source/PROCESSING/NOISEESTIMATION/sources.cmake)
include(source/PROCESSING/ID/sources.cmake)
include(source/PROCESSING/CENTROIDING/sources.cmake)
include(source/PROCESSING/SPECTRAMERGING/sources.cmake)
include(source/PROCESSING/SCALING/sources.cmake)
include(source/FEATUREFINDER/sources.cmake)
include(source/COMPARISON/sources.cmake)
include(source/QC/sources.cmake)
if(NOT DISABLE_OPENSWATH)
  include(source/ANALYSIS/OPENSWATH/sources.cmake)
  include(source/ANALYSIS/OPENSWATH/DATAACCESS/sources.cmake)
endif()
include(source/APPLICATIONS/sources.cmake)

set(OpenMS_Algo_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS Algo source files")
set(OpenMS_sources CACHE INTERNAL "")

## Header files
set(OpenMS_Algo_sources_h CACHE INTERNAL "OpenMS Algo header files")

include(include/OpenMS/ML/INTERPOLATION/sources.cmake)
include(include/OpenMS/ML/NNLS/sources.cmake)
include(include/OpenMS/ML/SVM/sources.cmake)
include(include/OpenMS/ML/RANSAC/sources.cmake)
include(include/OpenMS/ML/REGRESSION/sources.cmake)
include(include/OpenMS/ML/CLUSTERING/sources.cmake)
include(include/OpenMS/ML/GRIDSEARCH/sources.cmake)
include(include/OpenMS/ML/CROSSVALIDATION/sources.cmake)
include(include/OpenMS/ML/ROCCURVE/sources.cmake)
include(include/OpenMS/ANALYSIS/DECHARGING/sources.cmake)
include(include/OpenMS/ANALYSIS/ID/sources.cmake)
include(include/OpenMS/ANALYSIS/MAPMATCHING/sources.cmake)
include(include/OpenMS/ANALYSIS/QUANTITATION/sources.cmake)
include(include/OpenMS/ANALYSIS/SEQUENCE/sources.cmake)
include(include/OpenMS/ANALYSIS/MRM/sources.cmake)
include(include/OpenMS/ANALYSIS/NUXL/sources.cmake)
include(include/OpenMS/ANALYSIS/TARGETED/sources.cmake)
include(include/OpenMS/ANALYSIS/TOPDOWN/sources.cmake)
include(include/OpenMS/ANALYSIS/XLMS/sources.cmake)
include(include/OpenMS/COMPARISON/sources.cmake)
include(include/OpenMS/PROCESSING/BASELINE/sources.cmake)
include(include/OpenMS/PROCESSING/CALIBRATION/sources.cmake)
include(include/OpenMS/PROCESSING/MISC/sources.cmake)
include(include/OpenMS/PROCESSING/DEISOTOPING/sources.cmake)
include(include/OpenMS/PROCESSING/ID/sources.cmake)
include(include/OpenMS/PROCESSING/FEATURE/sources.cmake)
include(include/OpenMS/PROCESSING/NOISEESTIMATION/sources.cmake)
include(include/OpenMS/PROCESSING/SMOOTHING/sources.cmake)
include(include/OpenMS/PROCESSING/FILTERING/sources.cmake)
include(include/OpenMS/PROCESSING/RESAMPLING/sources.cmake)
include(include/OpenMS/PROCESSING/SCALING/sources.cmake)
include(include/OpenMS/FEATUREFINDER/sources.cmake)
include(include/OpenMS/PROCESSING/CENTROIDING/sources.cmake)
include(include/OpenMS/PROCESSING/SPECTRAMERGING/sources.cmake)
include(include/OpenMS/QC/sources.cmake)
if(NOT DISABLE_OPENSWATH)
  include(include/OpenMS/ANALYSIS/OPENSWATH/sources.cmake)
  include(include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/sources.cmake)
endif()
include(include/OpenMS/APPLICATIONS/sources.cmake)

set(OpenMS_Algo_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS Algo header files")
set(OpenMS_sources_h CACHE INTERNAL "")
```

**Step 2: Commit**

```bash
git add src/openms/algo_includes.cmake
git commit -m "build: add algo_includes.cmake for OpenMS_Algo library sources"
```

---

### Task 4: Handle DLL Export Macro (`OPENMS_DLLAPI`)

**Files:**
- Create: `src/openms/include/OpenMS/OpenMSDLLConfig.h`
- Modify: `src/openms/include/OpenMS/config.h.in` (line ~25)

**Context:** Currently, `OpenMSConfig.h` is auto-generated by `generate_export_header()` and defines `OPENMS_DLLAPI`. With three libraries, each would get its own export macro (`OPENMS_CORE_DLLAPI`, etc.). To avoid changing every header file, we create a shared export header that all three libraries use.

**Step 1: Create shared DLL export header**

```cpp
// OpenMSDLLConfig.h — Shared visibility/export macros for all OpenMS sub-libraries
//
// This header defines OPENMS_DLLAPI for use by all OpenMS sub-libraries.
// It replaces the auto-generated per-library export headers to maintain
// backward compatibility during the library decomposition.

#ifndef OPENMS_DLLCONFIG_H
#define OPENMS_DLLCONFIG_H

#ifdef _MSC_VER
  // Windows: export when building any OpenMS sub-library, import when consuming
  #if defined(OpenMS_Core_EXPORTS) || defined(OpenMS_IO_EXPORTS) || defined(OpenMS_Algo_EXPORTS) || defined(OpenMS_EXPORTS)
    #define OPENMS_DLLAPI __declspec(dllexport)
  #else
    #define OPENMS_DLLAPI __declspec(dllimport)
  #endif
  // Suppress DLL interface warnings for STL types
  #pragma warning(disable: 4251)
#else
  // Unix: use visibility attribute
  #if __GNUC__ >= 4
    #define OPENMS_DLLAPI __attribute__((visibility("default")))
  #else
    #define OPENMS_DLLAPI
  #endif
#endif

// Deprecated export macro — will be per-library in a future version
#ifndef OPENMS_DEPRECATED
  #ifdef _MSC_VER
    #define OPENMS_DEPRECATED __declspec(deprecated)
  #else
    #define OPENMS_DEPRECATED __attribute__((__deprecated__))
  #endif
#endif

#endif // OPENMS_DLLCONFIG_H
```

**Step 2: Update config.h.in to include the new header**

In `src/openms/include/OpenMS/config.h.in`, change line 25 from:
```cpp
#include <OpenMS/OpenMSConfig.h>
```
to:
```cpp
#include <OpenMS/OpenMSDLLConfig.h>
```

**Step 3: Modify openms_add_library to skip generate_export_header for sub-libraries**

In `cmake/add_library_macros.cmake`, the `DLL_EXPORT_PATH` parameter controls export header generation (lines 170-183). For sub-libraries, we pass an empty `DLL_EXPORT_PATH` to skip generation, since they use the shared header.

**Note:** The monolithic `OpenMS` target will no longer exist as a real library, so the old `OpenMSConfig.h` won't be generated. The new `OpenMSDLLConfig.h` replaces it.

**Step 4: Commit**

```bash
git add src/openms/include/OpenMS/OpenMSDLLConfig.h
git commit -m "build: add shared DLL export header for multi-library support"
```

---

### Task 5: Create sub-library CMakeLists.txt files

**Files:**
- Create: `src/openms/core/CMakeLists.txt`
- Create: `src/openms/io/CMakeLists.txt`
- Create: `src/openms/algo/CMakeLists.txt`
- Modify: `src/openms/CMakeLists.txt`

**Context:** Each sub-library gets its own CMakeLists.txt that calls `openms_add_library()`. The parent `src/openms/CMakeLists.txt` orchestrates all three and creates the umbrella target.

**Step 1: Create core/CMakeLists.txt**

```cmake
project("OpenMS_Core")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# Aggregate Core source files
include(${CMAKE_CURRENT_SOURCE_DIR}/../core_includes.cmake)

# Add configured headers
source_group("Header Files\\OpenMS" FILES ${OpenMS_configured_headers})
list(APPEND OpenMS_Core_sources ${OpenMS_Core_sources_h} ${OpenMS_configured_headers})
list(REMOVE_DUPLICATES OpenMS_Core_sources)

# External dependencies for Core
set(CORE_DEP_LIBRARIES
  LibSVM::LibSVM
  XercesC::XercesC
  Qt6::Core
)

set(CORE_DEP_PRIVATE_LIBRARIES
  ${LPTARGET}
  BZip2::BZip2
  Boost::boost
  Boost::regex
  Eigen3::Eigen
  GTE
  Quadtree
  SIMDe
  ZLIB::ZLIB
)

if (OPENMP_FOUND)
  list(APPEND CORE_DEP_LIBRARIES OpenMP::OpenMP_CXX)
endif()

if(APPLE)
  find_library(CoreFoundation_LIBRARY CoreFoundation)
  find_library(CoreServices_LIBRARY CoreServices)
  list(APPEND CORE_DEP_LIBRARIES ${CoreFoundation_LIBRARY} ${CoreServices_LIBRARY})
endif()

openms_add_library(TARGET_NAME  OpenMS_Core
                   SOURCE_FILES  ${OpenMS_Core_sources}
                   HEADER_FILES  ${OpenMS_Core_sources_h}
                                 ${OpenMS_configured_headers}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                     ${CMAKE_CURRENT_BINARY_DIR}/../include
                   LINK_LIBRARIES ${CORE_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${CORE_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "")

if (MSVC)
  target_compile_options(OpenMS_Core PRIVATE "/we4100" "/we4189")
endif()
```

**Step 2: Create io/CMakeLists.txt**

```cmake
project("OpenMS_IO")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

include(${CMAKE_CURRENT_SOURCE_DIR}/../io_includes.cmake)

source_group("Header Files\\OpenMS" FILES ${OpenMS_configured_headers})
list(APPEND OpenMS_IO_sources ${OpenMS_IO_sources_h})
list(REMOVE_DUPLICATES OpenMS_IO_sources)

set(IO_DEP_LIBRARIES
  OpenMS_Core
  Qt6::Network
)

set(IO_DEP_PRIVATE_LIBRARIES
  Boost::boost
  Boost::date_time
  Boost::iostreams
  Boost::regex
  Eigen3::Eigen
  SQLiteCpp
  nlohmann_json::nlohmann_json
  ZLIB::ZLIB
  BZip2::BZip2
  SIMDe
)

if (WITH_PARQUET)
  list(APPEND IO_DEP_PRIVATE_LIBRARIES
       ${OPENMS_ARROW_TARGET}
       ${OPENMS_ARROW_COMPUTE_TARGET}
       ${OPENMS_PARQUET_TARGET})
  if (OPENMS_ARROW_DATASET_TARGET)
    list(APPEND IO_DEP_PRIVATE_LIBRARIES ${OPENMS_ARROW_DATASET_TARGET} LibXml2::LibXml2)
  endif()
endif()

if (ENABLE_TDL)
  list(APPEND IO_DEP_PRIVATE_LIBRARIES tdl::tdl)
endif()

if (WITH_HDF5)
  list(APPEND IO_DEP_PRIVATE_LIBRARIES HDF5::HDF5)
endif()

openms_add_library(TARGET_NAME  OpenMS_IO
                   SOURCE_FILES  ${OpenMS_IO_sources}
                   HEADER_FILES  ${OpenMS_IO_sources_h}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                     ${CMAKE_CURRENT_BINARY_DIR}/../include
                   LINK_LIBRARIES ${IO_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${IO_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "")

if (ENABLE_TDL)
  target_compile_definitions(OpenMS_IO PUBLIC ENABLE_TDL)
endif()
if (WITH_PARQUET)
  target_compile_definitions(OpenMS_IO PUBLIC WITH_PARQUET=1)
endif()
if (OPENMS_ARROW_DATASET_TARGET)
  target_compile_definitions(OpenMS_IO PRIVATE WITH_ARROW_DATASET)
endif()

if (MSVC)
  target_compile_options(OpenMS_IO PRIVATE "/we4100" "/we4189")
endif()
```

**Step 3: Create algo/CMakeLists.txt**

```cmake
project("OpenMS_Algo")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

include(${CMAKE_CURRENT_SOURCE_DIR}/../algo_includes.cmake)

source_group("Header Files\\OpenMS" FILES ${OpenMS_configured_headers})
list(APPEND OpenMS_Algo_sources ${OpenMS_Algo_sources_h})
list(REMOVE_DUPLICATES OpenMS_Algo_sources)

set(ALGO_DEP_LIBRARIES
  OpenMS_IO
  OpenMS_Core
)

set(ALGO_DEP_PRIVATE_LIBRARIES
  OpenSwathAlgo
  Boost::boost
  Boost::regex
  Eigen3::Eigen
  Evergreen
  IsoSpec
  eol-bspline
)

openms_add_library(TARGET_NAME  OpenMS_Algo
                   SOURCE_FILES  ${OpenMS_Algo_sources}
                   HEADER_FILES  ${OpenMS_Algo_sources_h}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                     ${CMAKE_CURRENT_BINARY_DIR}/../include
                   LINK_LIBRARIES ${ALGO_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${ALGO_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "")

if (MSVC)
  target_compile_options(OpenMS_Algo PRIVATE "/we4100" "/we4189")
endif()
```

**Step 4: Modify parent src/openms/CMakeLists.txt**

Replace the monolithic library build with sub-library orchestration. Keep the extern subdirectory, configh.cmake config, and install commands. Replace the single `openms_add_library(OpenMS ...)` call with:

```cmake
# Build sub-libraries
add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(algo)

# Backward-compatible umbrella target
add_library(OpenMS INTERFACE)
target_link_libraries(OpenMS INTERFACE OpenMS_Algo OpenMS_IO OpenMS_Core)
target_include_directories(OpenMS INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
  "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>"
)
openms_register_export_target(OpenMS)
install(TARGETS OpenMS EXPORT OpenMSTargets)
```

Keep the `share/` and `cmake/Modules` install commands, doxygen path registration.

**Step 5: Build and validate**

```bash
cmake --build OpenMS-build -j$(nproc) 2>&1 | head -100
```

Expected: Three libraries build successfully. Watch for:
- Missing source files in any library
- Unresolved symbols (dependency in wrong library)
- DLL export issues

**Step 6: Run tests**

```bash
ctest --test-dir OpenMS-build --output-on-failure -j$(nproc)
```

Expected: All existing tests pass. Tests still link against `${OpenMS_LIBRARIES}` which will now resolve through the umbrella.

**Step 7: Commit**

```bash
git add src/openms/core/ src/openms/io/ src/openms/algo/ src/openms/CMakeLists.txt
git commit -m "build: split OpenMS into Core, IO, and Algo sub-libraries

Three libraries with strict layering: Algo → IO → Core.
Backward-compatible OpenMS INTERFACE target links all three.
No source file or header changes — only CMake wiring."
```

---

### Task 6: Fix compilation issues from library split

**Context:** The initial build will likely fail due to:
1. `sources.cmake` append pattern — each file appends to `OpenMS_sources`, which gets reset between includes.cmake files. Need to verify the reset-and-capture pattern works.
2. Missing dependencies — some modules may need libraries not assigned to their tier.
3. Symbol visibility — without `generate_export_header`, need to verify `OpenMSDLLConfig.h` is found.

**Step 1: Fix sources.cmake accumulator pattern**

The key issue: each `sources.cmake` does `set(OpenMS_sources ${OpenMS_sources} ${sources})`. After `core_includes.cmake` captures into `OpenMS_Core_sources` and resets `OpenMS_sources`, the next includes file (`io_includes.cmake`) starts fresh. This should work, but verify by adding debug messages:

```cmake
# At end of core_includes.cmake
message(STATUS "OpenMS_Core: ${OpenMS_Core_sources_count} source files")
list(LENGTH OpenMS_Core_sources OpenMS_Core_sources_count)
message(STATUS "OpenMS_Core source count: ${OpenMS_Core_sources_count}")
```

**Step 2: Fix any dependency misattributions**

Common issues to watch for:
- FORMAT files that need `Eigen3::Eigen` → already in IO's private deps
- ANALYSIS files that need `Qt6::Network` → add to Algo's deps if needed
- CHEMISTRY files that need `Boost::date_time` → add to Core's private deps if needed

For each linker error, identify which library the missing symbol lives in and add the appropriate `target_link_libraries` entry.

**Step 3: Fix OPENMS_DLLAPI resolution**

Verify that `config.h` → `OpenMSDLLConfig.h` include chain works. Each sub-library's compilation should see `OPENMS_DLLAPI` defined. Check with:
```bash
# Should show the define is present
grep -r "OpenMS_Core_EXPORTS\|OpenMS_IO_EXPORTS\|OpenMS_Algo_EXPORTS" OpenMS-build/
```

CMake automatically defines `<target>_EXPORTS` when building a shared library target. Verify these appear in compile flags.

**Step 4: Iterate until clean build**

Repeat build-fix-build until `cmake --build OpenMS-build -j$(nproc)` succeeds with no errors.

**Step 5: Run full test suite**

```bash
ctest --test-dir OpenMS-build --output-on-failure -j$(nproc)
```

**Step 6: Commit all fixes**

```bash
git add -A
git commit -m "build: fix compilation issues from library split"
```

---

### Task 7: Update test CMakeLists.txt for sub-library linking

**Files:**
- Modify: `src/tests/class_tests/openms/CMakeLists.txt`

**Context:** Tests currently link against `${OpenMS_LIBRARIES}` set by `openms_add_library(OpenMS ...)`. With the umbrella INTERFACE target, we need to update tests to link against `OpenMS` (the umbrella) or the specific sub-library.

**Step 1: Update test linking to use umbrella target**

In `src/tests/class_tests/openms/CMakeLists.txt`, line 42:
```cmake
# Old:
target_link_libraries(${_class_test} ${OpenMS_LIBRARIES})
# New:
target_link_libraries(${_class_test} OpenMS)
```

Also update `OpenMS_INCLUDE_DIRECTORIES` usage (line 27):
```cmake
# Old:
set(OPENMS_CLASS_TESTS_EXTERNAL_INCLUDE_DIRECTORIES "${OpenMS_INCLUDE_DIRECTORIES};${Boost_INCLUDE_DIRS}")
# New — includes propagate through target_link_libraries, so this may no longer be needed.
# But keep Boost for tests that directly use Boost headers:
set(OPENMS_CLASS_TESTS_EXTERNAL_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
```

**Step 2: Build and run tests**

```bash
cmake --build OpenMS-build -j$(nproc)
ctest --test-dir OpenMS-build --output-on-failure -j$(nproc)
```

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/CMakeLists.txt
git commit -m "test: update class tests to link against OpenMS umbrella target"
```

---

## Phase 2: Boundary Fixes

### Task 8: Move FASTAContainer from DATASTRUCTURES to FORMAT

**Files:**
- Move: `src/openms/source/DATASTRUCTURES/FASTAContainer.cpp` → `src/openms/source/FORMAT/FASTAContainer.cpp` (if it exists)
- Move: `src/openms/include/OpenMS/DATASTRUCTURES/FASTAContainer.h` → `src/openms/include/OpenMS/FORMAT/FASTAContainer.h`
- Modify: `src/openms/source/DATASTRUCTURES/sources.cmake` (remove FASTAContainer)
- Modify: `src/openms/source/FORMAT/sources.cmake` (add FASTAContainer)
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/sources.cmake` (remove)
- Modify: `src/openms/include/OpenMS/FORMAT/sources.cmake` (add)
- Add: `src/openms/include/OpenMS/DATASTRUCTURES/FASTAContainer.h` — compatibility redirect header

**Step 1: Move the files**

```bash
git mv src/openms/include/OpenMS/DATASTRUCTURES/FASTAContainer.h src/openms/include/OpenMS/FORMAT/FASTAContainer.h
# Check if .cpp exists:
ls src/openms/source/DATASTRUCTURES/FASTAContainer.cpp
# If it exists, move it too
```

**Step 2: Create compatibility redirect header**

At the old location, create a header that includes the new location:
```cpp
// Compatibility redirect — FASTAContainer moved to FORMAT module
#pragma once
#include <OpenMS/FORMAT/FASTAContainer.h>
```

**Step 3: Update sources.cmake files**

Remove from DATASTRUCTURES `sources.cmake`, add to FORMAT `sources.cmake`.

**Step 4: Update any #include directives in the codebase**

```bash
grep -r "DATASTRUCTURES/FASTAContainer" src/openms/ --include="*.cpp" --include="*.h"
```

Update these to use `FORMAT/FASTAContainer.h`.

**Step 5: Build and test**

```bash
cmake --build OpenMS-build -j$(nproc)
ctest --test-dir OpenMS-build --output-on-failure -j$(nproc)
```

**Step 6: Commit**

```bash
git add -A
git commit -m "refactor: move FASTAContainer from DATASTRUCTURES to FORMAT

FASTAContainer depends on FASTAFile and is format-specific I/O.
Moving it to FORMAT places it in the IO tier of the library split.
A compatibility redirect header is left at the old location."
```

---

### Task 9: Move QTCluster from DATASTRUCTURES to Algo tier

**Files:**
- Move: `src/openms/include/OpenMS/DATASTRUCTURES/QTCluster.h` → `src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/QTCluster.h` (or FEATUREFINDER)
- Move: `src/openms/source/DATASTRUCTURES/QTCluster.cpp` → matching source location
- Update: sources.cmake files
- Add: compatibility redirect at old location
- Update: all #include references

**Context:** QTCluster depends on `CHEMISTRY/AASequence.h`, which would be fine in Core since CHEMISTRY is in Core. However, QTCluster is algorithmically a clustering data structure used by QTClusterFinder in ANALYSIS/MAPMATCHING. Its dependency on AASequence is domain-specific (peptide clustering), making it an algorithm component.

**Step 1: Move files**

```bash
git mv src/openms/include/OpenMS/DATASTRUCTURES/QTCluster.h src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/QTCluster.h
git mv src/openms/source/DATASTRUCTURES/QTCluster.cpp src/openms/source/ANALYSIS/MAPMATCHING/QTCluster.cpp
```

**Step 2: Create compatibility redirect**

**Step 3: Update sources.cmake files**

**Step 4: Update #include references**

```bash
grep -r "DATASTRUCTURES/QTCluster" src/ --include="*.cpp" --include="*.h"
```

**Step 5: Build, test, commit**

```bash
cmake --build OpenMS-build -j$(nproc)
ctest --test-dir OpenMS-build --output-on-failure -j$(nproc)
git add -A
git commit -m "refactor: move QTCluster from DATASTRUCTURES to ANALYSIS/MAPMATCHING

QTCluster depends on AASequence and is used only by QTClusterFinder.
Moving it to ANALYSIS/MAPMATCHING places it in the Algo tier."
```

---

### Task 10: Move PosteriorErrorProbabilityModel from MATH to Algo

**Files:**
- Move: `src/openms/include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h` → `src/openms/include/OpenMS/ANALYSIS/ID/PosteriorErrorProbabilityModel.h`
- Move: `src/openms/source/MATH/STATISTICS/PosteriorErrorProbabilityModel.cpp` → `src/openms/source/ANALYSIS/ID/PosteriorErrorProbabilityModel.cpp`
- Update: sources.cmake files
- Add: compatibility redirect
- Update: all #include references

**Context:** PosteriorErrorProbabilityModel depends on METADATA peptide identification types. It's a statistical model for peptide scoring, which is an analysis algorithm.

Follow same pattern as Tasks 8-9.

---

### Task 11: Refactor OnDiscMSExperiment (Core/IO boundary)

**Files:**
- Modify: `src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h`
- Modify: `src/openms/source/KERNEL/OnDiscMSExperiment.cpp`
- Possibly create: abstract base or pimpl implementation

**Context:** OnDiscMSExperiment is in KERNEL (Core) but depends on `FORMAT/HANDLERS/IndexedMzMLHandler.h` and `FORMAT/OPTIONS/PeakFileOptions.h` (IO). This creates a Core→IO back-dependency.

**Approach options (decide during implementation):**
1. **Pimpl pattern**: Move the FORMAT-dependent members behind an opaque pointer. Core header has no FORMAT includes. IO library provides the implementation.
2. **Move to IO**: If OnDiscMSExperiment is primarily used by IO code (FileHandler, IndexedMzMLFileLoader), it could live in IO with a typedef/alias in KERNEL for compatibility.
3. **Forward declarations**: Replace `#include` with forward declarations in the header, move includes to .cpp (but .cpp is in Core, which still can't see IO).

**Recommended approach (pimpl):**

```cpp
// OnDiscMSExperiment.h (in KERNEL/Core)
// Remove these includes:
// #include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>
// #include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

// Add forward declaration:
namespace OpenMS { namespace Internal { class OnDiscMSExperimentImpl; } }

class OPENMS_DLLAPI OnDiscMSExperiment
{
public:
    // ... existing public interface ...
private:
    std::unique_ptr<Internal::OnDiscMSExperimentImpl> impl_;
};
```

The implementation class `OnDiscMSExperimentImpl` lives in IO and is registered/created through a factory or constructor that takes a filename.

**Note:** This is the most complex refactoring task. It may require careful API design to maintain backward compatibility. The exact approach should be decided based on how OnDiscMSExperiment is used across the codebase.

**Step 1: Analyze usage patterns**

```bash
grep -r "OnDiscMSExperiment" src/ --include="*.cpp" --include="*.h" -l
```

**Step 2: Implement chosen refactoring approach**

**Step 3: Build and test thoroughly**

**Step 4: Commit**

---

### Task 12: Move specialized FORMAT handlers to Algo tier

**Files to move from FORMAT to appropriate ANALYSIS subdirectories:**

| File | Destination | Reason |
|------|-------------|--------|
| `TraMLHandler.h/.cpp` | `ANALYSIS/TARGETED/` | Serializes TargetedExperiment |
| `TraMLFile.h/.cpp` | `ANALYSIS/TARGETED/` | Serializes TargetedExperiment |
| `MRMFeatureQCFile.h/.cpp` | `ANALYSIS/OPENSWATH/` | OpenSWATH-specific |
| `MRMFeaturePickerFile.h/.cpp` | `ANALYSIS/OPENSWATH/` | OpenSWATH-specific |
| `XQuestResultXMLFile.h/.cpp` | `ANALYSIS/XLMS/` | Cross-linking specific |
| `FLASHDeconvFeatureFile.h/.cpp` | `ANALYSIS/TOPDOWN/` | TopDown specific |
| `FLASHDeconvSpectrumFile.h/.cpp` | `ANALYSIS/TOPDOWN/` | TopDown specific |
| `SwathFileConsumer.h/.cpp` | `ANALYSIS/OPENSWATH/` | SWATH-specific |
| `AbsoluteQuantitationMethodFile.h/.cpp` | `ANALYSIS/QUANTITATION/` | Quantitation-specific |
| `MRMFile.h/.cpp` | `ANALYSIS/MRM/` | MRM-specific |
| `IBSpectraFile.h/.cpp` | `ANALYSIS/QUANTITATION/` | Isobaric-specific |
| `GNPSMGFFile.h/.cpp` | `ANALYSIS/ID/` | Depends on SpectraMerger |

**For each file:**
1. `git mv` source and header to new location
2. Create compatibility redirect at old location
3. Update sources.cmake (remove from FORMAT, add to ANALYSIS/*)
4. Update #include references
5. Build and test

**Tip:** Do this in batches by destination directory (all OPENSWATH files together, all TARGETED together, etc.).

**After all moves: Build and test**

```bash
cmake --build OpenMS-build -j$(nproc)
ctest --test-dir OpenMS-build --output-on-failure -j$(nproc)
```

**Commit after each batch or after all moves.**

---

## Phase 3: Test System Split

### Task 13: Create per-library test directories and CMakeLists

**Files:**
- Create: `src/tests/class_tests/openms_core/CMakeLists.txt`
- Create: `src/tests/class_tests/openms_core/executables.cmake`
- Create: `src/tests/class_tests/openms_io/CMakeLists.txt`
- Create: `src/tests/class_tests/openms_io/executables.cmake`
- Create: `src/tests/class_tests/openms_algo/CMakeLists.txt`
- Create: `src/tests/class_tests/openms_algo/executables.cmake`
- Modify: `src/tests/CMakeLists.txt` (add new subdirectories)

**Context:** Current test structure has all tests in `src/tests/class_tests/openms/` with `executables.cmake` defining per-module test lists. We split these lists into per-library executables.cmake files.

**Step 1: Create openms_core test infrastructure**

Core test executables (from existing `executables.cmake`):
```cmake
# openms_core/executables.cmake
set(TEST_executables
    ${concept_executables_list}
    ${system_executables_list}
    ${datastructures_executables_list}
    ${kernel_executables_list}
    ${metadata_executables_list}
    ${chemistry_executables_list}
    ${math_executables_list}
    ${ionmobility_executables_list}
)
```

Core `CMakeLists.txt` — similar to existing but links `OpenMS_Core` only:
```cmake
# Each test links only against OpenMS_Core
target_link_libraries(${_class_test} OpenMS_Core)
```

**Step 2: Create openms_io and openms_algo test infrastructure**

IO tests:
```cmake
set(TEST_executables ${format_executables_list})
# Link: OpenMS_IO
```

Algo tests:
```cmake
set(TEST_executables
    ${filtering_executables_list}
    ${comparison_executables_list}
    ${analysis_executables_list}
    ${applications_executables_list}
    ${transformations_executables_list}
    ${swath_executables_list}
    ${qc_executables_list}
)
# Link: OpenMS_Algo
```

**Step 3: Add CTest labels**

In each new CMakeLists.txt, after `add_test()`:
```cmake
set_tests_properties(${_class_test} PROPERTIES LABELS "core")  # or "io", "algo"
```

**Step 4: Symlink or reference shared test data**

The test data directory `src/tests/class_tests/openms/data/` should be accessible to all three test dirs. Use `CF_OPENMS_TEST_DATA_PATH` pointing to the shared data directory.

**Step 5: Build and test**

```bash
cmake --build OpenMS-build -j$(nproc)
ctest --test-dir OpenMS-build -L core --output-on-failure
ctest --test-dir OpenMS-build -L io --output-on-failure
ctest --test-dir OpenMS-build -L algo --output-on-failure
```

**Step 6: Commit**

```bash
git add -A
git commit -m "test: split class tests into per-library test directories

Tests now organized by library: openms_core, openms_io, openms_algo.
Each test dir links only against its own library, enforcing boundaries.
CTest labels allow running per-library: ctest -L core"
```

---

### Task 14: Fix test compilation issues

**Context:** Core tests that only link `OpenMS_Core` may fail if they indirectly depend on IO or Algo functionality. These test failures are valuable — they reveal actual boundary violations.

**For each failing test:**
1. Check what symbols are missing
2. Determine if the test genuinely belongs in a different tier
3. Move the test to the correct tier, or fix the dependency

**Common patterns:**
- `OnDiscMSExperiment_test` may need IO → move to IO tests (until pimpl refactoring is done)
- `FASTAContainer_test` → already in IO after Task 8
- Tests for moved files (QTCluster_test, PosteriorErrorProbabilityModel_test) → move to Algo tests

---

## Phase 4: Cleanup and Validation

### Task 15: Update OpenMSConfig.cmake.in for multi-library export

**Files:**
- Modify: `cmake/OpenMSConfig.cmake.in`

**Context:** External consumers need to find all three libraries plus the umbrella. The config file template needs updating.

**Step 1: Update config template**

Add dependency finding for all three libraries and the umbrella target. Since the umbrella is INTERFACE and links all three, consumers just need `find_package(OpenMS)` and `target_link_libraries(... OpenMS)`.

Ensure `OpenMSTargets.cmake` includes all four targets (Core, IO, Algo, OpenMS).

**Step 2: Test external consumption**

Create a minimal CMake project that does:
```cmake
find_package(OpenMS REQUIRED)
add_executable(test_consumer main.cpp)
target_link_libraries(test_consumer OpenMS)
```

---

### Task 16: Update pyOpenMS CMakeLists.txt

**Files:**
- Modify: `src/pyOpenMS/CMakeLists.txt`

**Context:** pyOpenMS links against the OpenMS library. Update to use the umbrella target.

**Step 1: Update linking**

Replace any reference to `${OpenMS_LIBRARIES}` or specific OpenMS target with the `OpenMS` umbrella target.

**Step 2: Build and test pyOpenMS**

```bash
cmake --build OpenMS-build --target pyopenms -j$(nproc)
PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/ -v
```

---

### Task 17: Full validation

**Step 1: Clean build from scratch**

```bash
rm -rf OpenMS-build && mkdir OpenMS-build && cd OpenMS-build
cmake .. [your usual cmake flags]
cmake --build . -j$(nproc)
```

**Step 2: Run all tests**

```bash
ctest --test-dir OpenMS-build --output-on-failure -j$(nproc)
```

**Step 3: Run per-library tests to verify isolation**

```bash
ctest --test-dir OpenMS-build -L core --output-on-failure
ctest --test-dir OpenMS-build -L io --output-on-failure
ctest --test-dir OpenMS-build -L algo --output-on-failure
```

**Step 4: Build pyOpenMS and run Python tests**

```bash
cmake --build OpenMS-build --target pyopenms -j$(nproc)
PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/ -v
```

**Step 5: Verify no regressions**

Compare test results against the pre-split baseline. All tests that passed before should still pass.

---

### Task 18: Remove legacy monolithic build path

**Files:**
- Remove: `src/openms/includes.cmake` (replaced by core/io/algo_includes.cmake)
- Clean up: `src/tests/class_tests/openms/` (if all tests moved to per-library dirs)

**Only do this after all tests pass and the split is validated.**

---

## Known Issues and Future Work

1. **Windows DLL export**: The shared `OpenMSDLLConfig.h` uses a simplified approach where all sub-library builds define the same `OPENMS_DLLAPI` as `dllexport`. This works but isn't technically correct for cross-library header inclusion. Future work: per-library export macros (`OPENMS_CORE_DLLAPI`, `OPENMS_IO_DLLAPI`, `OPENMS_ALGO_DLLAPI`).

2. **Specialized format handler tests**: Tests for moved format handlers (TraMLFile_test, etc.) need to move to the algo test directory.

3. **TOPP tools**: Currently link against monolithic OpenMS. Update to use umbrella target (trivial change, same behavior).

4. **Unity build**: May need updating to work with three separate libraries.

5. **CI pipeline**: Add per-library build/test stages for faster CI feedback.
