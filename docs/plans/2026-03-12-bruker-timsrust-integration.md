# Bruker TimsTOF Integration Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add native Bruker TimsTOF `.d` file reading via timsrust_cpp_bridge, supporting DDA-PASEF, DIA-PASEF, and raw frame-level 4D access.

**Architecture:** A new `BrukerTimsFile` class wraps the timsrust C ABI to load `.d` directories into `MSExperiment`. It's registered as `BRUKER_TDF` in `FileTypes`, dispatched by `FileHandler`, and exposed through `FileConverter` with timsrust-specific parameters. The entire integration is gated behind `WITH_TIMSRUST` (CMake, default ON).

**Tech Stack:** C++, CMake FetchContent, timsrust_cpp_bridge (C ABI static library), OpenMS FORMAT/KERNEL/IONMOBILITY modules.

**Spec:** `docs/specs/2026-03-12-bruker-timsrust-integration-design.md`

---

## Chunk 1: CMake and File Type Infrastructure

### Task 1: CMake Option and FetchContent

**Files:**
- Modify: `CMakeLists.txt:60` (after `WITH_PARQUET` option)
- Modify: `cmake/cmake_findExternalLibs.cmake:178` (after Arrow/Parquet section)
- Modify: `src/openms/CMakeLists.txt:89-100` (after Parquet linking block)
- Modify: `src/openms/CMakeLists.txt:164-167` (after Parquet compile definition)
- Modify: `cmake/OpenMSConfig.cmake.in:72` (after `WITH_PARQUET`)

- [ ] **Step 1: Add CMake options**

In `CMakeLists.txt`, after line 60 (`WITH_PARQUET`):

```cmake
option(WITH_TIMSRUST "Build with Bruker TimsTOF .d file support via timsrust" ON)
option(ENABLE_TIMSRUST_TESTS "Download Bruker test data for timsrust integration tests" OFF)
```

- [ ] **Step 2: Add FetchContent in cmake_findExternalLibs.cmake**

After the Arrow/Parquet `endif()` (after line 262 — must be outside the `if (WITH_PARQUET)` block), add:

```cmake
#------------------------------------------------------------------------------
# timsrust_cpp_bridge (Bruker TimsTOF .d file reading)
if (WITH_TIMSRUST)
  include(FetchContent)
  set(TIMSRUST_VERSION "0.1.0" CACHE STRING "timsrust_cpp_bridge version to fetch")

  # Detect platform for archive selection
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|ARM64")
      set(_TIMSRUST_PLATFORM "linux-aarch64")
    else()
      set(_TIMSRUST_PLATFORM "linux-x86_64")
    endif()
  elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
      set(_TIMSRUST_PLATFORM "macos-arm64")
    else()
      message(WARNING "timsrust_cpp_bridge: no pre-built Intel macOS binaries; build from source or use arm64")
      set(_TIMSRUST_PLATFORM "macos-arm64")  # fallback, may fail at link time
    endif()
  elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(_TIMSRUST_PLATFORM "windows-x86_64")
  else()
    message(FATAL_ERROR "Unsupported platform for timsrust_cpp_bridge: ${CMAKE_SYSTEM_NAME}")
  endif()

  set(_TIMSRUST_URL "https://github.com/OpenMS/timsrust_cpp_bridge/releases/download/v${TIMSRUST_VERSION}/timsrust_cpp_bridge-${_TIMSRUST_PLATFORM}.tar.gz")

  FetchContent_Declare(
    timsrust_cpp_bridge
    URL ${_TIMSRUST_URL}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  # Use Populate (not MakeAvailable) since the archive is pre-built with no CMakeLists.txt
  FetchContent_GetProperties(timsrust_cpp_bridge)
  if(NOT timsrust_cpp_bridge_POPULATED)
    FetchContent_Populate(timsrust_cpp_bridge)
  endif()

  # The extracted archive contains lib/cmake/ with config files
  list(APPEND CMAKE_PREFIX_PATH "${timsrust_cpp_bridge_SOURCE_DIR}")
  find_package(timsrust_cpp_bridge REQUIRED)
  message(STATUS "Found timsrust_cpp_bridge: ${timsrust_cpp_bridge_SOURCE_DIR}")
endif()
```

- [ ] **Step 3: Add conditional linking in src/openms/CMakeLists.txt**

Note: We use `target_compile_definitions(OpenMS PUBLIC WITH_TIMSRUST=1)` (Step 4) instead of `#cmakedefine` in `config.h.in`, matching the existing `WITH_PARQUET` pattern. No `config.h.in` change needed.

After the `WITH_PARQUET` linking block (after line 100):

```cmake
if (WITH_TIMSRUST)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES timsrust_cpp_bridge::timsrust_cpp_bridge)
endif()
```

After the `WITH_PARQUET` compile definition (after line 167):

```cmake
if (WITH_TIMSRUST)
    target_compile_definitions(OpenMS PUBLIC WITH_TIMSRUST=1)
endif()
```

- [ ] **Step 4: Export in OpenMSConfig.cmake.in**

In `cmake/OpenMSConfig.cmake.in`, after line 72 (`WITH_PARQUET`):

```cmake
set(WITH_TIMSRUST @WITH_TIMSRUST@)
```

- [ ] **Step 5: Verify CMake configuration**

Run: `cmake -B OpenMS-build -DWITH_TIMSRUST=OFF` (OFF for now since no release exists yet)
Expected: Configuration succeeds, `WITH_TIMSRUST` option visible in CMake output.

- [ ] **Step 6: Commit**

```bash
git add CMakeLists.txt cmake/cmake_findExternalLibs.cmake cmake/OpenMSConfig.cmake.in \
  src/openms/CMakeLists.txt
git commit -m "build: add WITH_TIMSRUST CMake option and FetchContent for timsrust_cpp_bridge"
```

---

### Task 2: File Type Registration (BRUKER_TDF)

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/FileTypes.h:95` (before `SIZE_OF_TYPE`)
- Modify: `src/openms/source/FORMAT/FileTypes.cpp:107` (before XML entry in array)
- Modify: `src/openms/source/FORMAT/FileTypes.cpp:251` (in `typeToMZML` switch)
- Modify: `src/tests/class_tests/openms/source/FileTypes_test.cpp`

- [ ] **Step 1: Write failing test for BRUKER_TDF type**

In `src/tests/class_tests/openms/source/FileTypes_test.cpp`, add after the existing `nameToType` tests for PARQUET (search for `PARQUET`):

```cpp
TEST_EQUAL(FileTypes::typeToName(FileTypes::BRUKER_TDF), "d");
TEST_EQUAL(FileTypes::BRUKER_TDF, FileTypes::nameToType("d"));
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cmake --build OpenMS-build --target FileTypes_test -j$(nproc) 2>&1 | tail -5`
Expected: Compilation error — `BRUKER_TDF` is not a member of `FileTypes`.

- [ ] **Step 3: Add BRUKER_TDF to FileTypes.h enum**

In `src/openms/include/OpenMS/FORMAT/FileTypes.h`, insert before `SIZE_OF_TYPE` (line 96):

```cpp
      BRUKER_TDF,         ///< Bruker TimsTOF .d directory (TDF format)
```

- [ ] **Step 4: Add TypeNameBinding in FileTypes.cpp**

In `src/openms/source/FORMAT/FileTypes.cpp`, insert before the XML entry (line 108, the comment says "make sure this comes last"):

```cpp
    TypeNameBinding(FileTypes::BRUKER_TDF, "d", "Bruker TDF", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE}),
```

- [ ] **Step 5: Add typeToMZML mapping**

In `src/openms/source/FORMAT/FileTypes.cpp`, in the `typeToMZML` switch (before `default:` at line 252):

```cpp
      case FileTypes::BRUKER_TDF: return "Bruker TDF format";
```

- [ ] **Step 6: Update FileTypes_test.cpp counts**

Find the `TEST_EQUAL` lines that check type counts (around line 133-135). Increment each by 1:
- `TEST_EQUAL(g.getTypes().size(), 41)` → change `41` to `42`
- `TEST_EQUAL(FileTypeList::typesWithProperties({}).size(), 64)` → change `64` to `65`

Also update the `SIZE_OF_TYPE` test if one exists.

- [ ] **Step 7: Build and run test**

Run:
```bash
cmake --build OpenMS-build --target FileTypes_test -j$(nproc)
ctest --test-dir OpenMS-build -R FileTypes_test -V
```
Expected: All FileTypes tests pass, including the new `BRUKER_TDF` assertions.

- [ ] **Step 8: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/FileTypes.h \
  src/openms/source/FORMAT/FileTypes.cpp \
  src/tests/class_tests/openms/source/FileTypes_test.cpp
git commit -m "feat: register BRUKER_TDF file type for Bruker TimsTOF .d directories"
```

---

### Task 3: FileHandler Directory Detection and Dispatch Stub

**Files:**
- Modify: `src/openms/source/FORMAT/FileHandler.cpp:161-169` (`getType()`)
- Modify: `src/openms/source/FORMAT/FileHandler.cpp:873` (`loadExperiment()` switch)
- Modify: `src/openms/source/FORMAT/FileHandler.cpp:895-907` (source file metadata)

- [ ] **Step 1: Add directory short-circuit in getType()**

In `src/openms/source/FORMAT/FileHandler.cpp`, modify `getType()` (lines 161-169). The key change: validate BRUKER_TDF **before** falling through to `getTypeByContent()`, which would fail on a directory path:

```cpp
  FileTypes::Type FileHandler::getType(const String& filename)
  {
    FileTypes::Type type = getTypeByFileName(filename);

    // Directory-based formats: validate marker files before returning.
    // Must happen before getTypeByContent() which would fail on directories.
    if (type == FileTypes::BRUKER_TDF)
    {
      String path = filename;
      // Strip trailing separator if present
      while (path.hasSuffix("/") || path.hasSuffix("\\"))
      {
        path = path.prefix(path.size() - 1);
      }
      // Validate directory contains TDF marker files
      if (File::exists(path + "/analysis.tdf") || File::exists(path + "/analysis.tdf_bin"))
      {
        return FileTypes::BRUKER_TDF;
      }
      return FileTypes::UNKNOWN; // .d suffix but not a TDF directory
    }

    if (type == FileTypes::UNKNOWN)
    {
      type = getTypeByContent(filename);
    }
    return type;
  }
```

- [ ] **Step 2: Add loadExperiment dispatch stub (compile-guarded)**

In `src/openms/source/FORMAT/FileHandler.cpp`, in the `loadExperiment()` switch, before the `default:` case (line 874), add:

```cpp
#ifdef WITH_TIMSRUST
      case FileTypes::BRUKER_TDF:
      {
        // TODO: BrukerTimsFile().load(filename, exp);
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      break;
#endif
```

- [ ] **Step 3: Handle directory paths in source file metadata**

In the `rewrite_source_file` block (lines 879-912), the code calls `File::basename(filename)` (line 895) and `computeFileHash(filename)` (line 907). For directories, `basename` may return empty string with trailing slashes, and `computeFileHash` would fail. Add a guard before the hash computation:

```cpp
      if (compute_hash && type != FileTypes::BRUKER_TDF)
      {
        src_file.setChecksum(computeFileHash(filename), SourceFile::ChecksumType::SHA1);
      }
```

And normalize the filename for `basename` and `setPathToFile`. Replace the existing `src_file.setNameOfFile(File::basename(filename))` (line 895) with:

```cpp
      // Normalize directory paths: strip trailing slashes for basename and URI
      String normalized_name = filename;
      while (normalized_name.hasSuffix("/") || normalized_name.hasSuffix("\\"))
        normalized_name = normalized_name.prefix(normalized_name.size() - 1);
      src_file.setNameOfFile(File::basename(normalized_name));
      String path_to_file = File::path(File::absolutePath(normalized_name));
```

- [ ] **Step 4: Build and verify**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: Builds successfully (the BRUKER_TDF case is behind `#ifdef WITH_TIMSRUST` which is off for now).

- [ ] **Step 5: Commit**

```bash
git add src/openms/source/FORMAT/FileHandler.cpp
git commit -m "feat: add directory-aware BRUKER_TDF detection and dispatch in FileHandler"
```

---

## Chunk 2: BrukerTimsFile Core Reader

### Task 4: BrukerTimsFile Header and Source Skeleton

**Files:**
- Create: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`
- Create: `src/openms/source/FORMAT/BrukerTimsFile.cpp`
- Modify: `src/openms/source/FORMAT/sources.cmake:121` (after WITH_PARQUET block)
- Modify: `src/openms/include/OpenMS/FORMAT/sources.cmake:133` (after WITH_PARQUET block)

- [ ] **Step 1: Create header file**

Create `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <cstdint>

namespace OpenMS
{

  /**
   * @brief Reader for Bruker TimsTOF .d directories via timsrust_cpp_bridge.
   *
   * Supports DDA-PASEF, DIA-PASEF, and raw frame-level 4D access.
   * Ion mobility data is stored in VSSC (1/K0) units using CONCATENATED format
   * for MS1 and DIA MS2, and scalar drift times for DDA MS2.
   */
  class OPENMS_DLLAPI BrukerTimsFile : public ProgressLogger
  {
  public:
    /// Processing and export configuration
    struct Config
    {
      uint32_t smoothing_window = 0;       ///< 0 = timsrust default; bin count
      uint32_t centroiding_window = 0;     ///< 0 = timsrust default; bin count
      double calibration_tolerance = 0.0;  ///< 0 = timsrust default; m/z tolerance
      bool calibrate = false;              ///< off by default (timsrust 0.4.2 caveat)

      enum ExportMode { AUTO, SPECTRUM, FRAME };
      ExportMode export_mode = AUTO;       ///< AUTO detects DDA vs DIA
    };

    /// Load entire .d directory into MSExperiment
    void load(const String& path, MSExperiment& exp, const Config& config = {});

    /// Streaming: read .d and feed spectra to consumer
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config = {});

  private:
    /// Load DDA-PASEF data: MS1 frames (CONCATENATED) + MS2 spectra (scalar IM)
    void loadDDA_(void* handle, MSExperiment& exp);

    /// Load DIA-PASEF data: MS1 frames (CONCATENATED) + MS2 frames split by SWATH window
    void loadDIA_(void* handle, MSExperiment& exp);

    /// Load raw frames without precursor grouping or SWATH splitting
    void loadFrames_(void* handle, MSExperiment& exp);

    /// Convert a single tims_frame to an MSSpectrum in CONCATENATED format
    void frameToSpectrum_(void* handle, const void* frame, MSSpectrum& spec) const;

    /// Detect DDA vs DIA by checking for SWATH windows
    bool isDIA_(void* handle) const;
  };

} // namespace OpenMS

#endif // WITH_TIMSRUST
```

- [ ] **Step 2: Create source skeleton**

Create `src/openms/source/FORMAT/BrukerTimsFile.cpp`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <timsrust_cpp_bridge.h>

#include <memory>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

namespace OpenMS
{

  // RAII wrapper for tims_dataset handle
  struct TimsDatasetHandle
  {
    tims_dataset* ptr = nullptr;

    ~TimsDatasetHandle()
    {
      if (ptr) tims_close(ptr);
    }

    TimsDatasetHandle() = default;
    TimsDatasetHandle(const TimsDatasetHandle&) = delete;
    TimsDatasetHandle& operator=(const TimsDatasetHandle&) = delete;
  };

  // RAII wrapper for tims_config
  struct TimsConfigHandle
  {
    tims_config* ptr = nullptr;

    ~TimsConfigHandle()
    {
      if (ptr) tims_config_free(ptr);
    }

    TimsConfigHandle() = default;
    TimsConfigHandle(const TimsConfigHandle&) = delete;
    TimsConfigHandle& operator=(const TimsConfigHandle&) = delete;
  };

  // Helper: get error string from handle
  static String getTimsError(tims_dataset* handle)
  {
    char buf[512];
    tims_get_last_error(handle, buf, sizeof(buf));
    return String(buf);
  }

  // Helper: open dataset with optional config
  static void openDataset(const String& path, const BrukerTimsFile::Config& config,
                          TimsDatasetHandle& ds, TimsConfigHandle& cfg)
  {
    bool has_config = (config.smoothing_window != 0 || config.centroiding_window != 0 ||
                       config.calibration_tolerance != 0.0 || config.calibrate);

    if (has_config)
    {
      cfg.ptr = tims_config_create();
      if (!cfg.ptr)
        throw Exception::BaseException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Internal error", "Failed to create timsrust config");

      if (config.smoothing_window != 0)
        tims_config_set_smoothing_window(cfg.ptr, config.smoothing_window);
      if (config.centroiding_window != 0)
        tims_config_set_centroiding_window(cfg.ptr, config.centroiding_window);
      if (config.calibration_tolerance != 0.0)
        tims_config_set_calibration_tolerance(cfg.ptr, config.calibration_tolerance);
      if (config.calibrate)
        tims_config_set_calibrate(cfg.ptr, true);

      timsffi_status status = tims_open_with_config(path.c_str(), cfg.ptr, &ds.ptr);
      if (status != TIMSFFI_OK)
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path,
          "timsrust error: " + getTimsError(ds.ptr));
    }
    else
    {
      timsffi_status status = tims_open(path.c_str(), &ds.ptr);
      if (status != TIMSFFI_OK)
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path,
          "timsrust error: " + getTimsError(ds.ptr));
    }
  }

  bool BrukerTimsFile::isDIA_(void* handle) const
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);
    unsigned int count = 0;
    tims_swath_window* windows = nullptr;
    timsffi_status status = tims_get_swath_windows(ds, &count, &windows);
    if (status == TIMSFFI_OK && windows)
    {
      tims_free_swath_windows(ds, windows);
    }
    return (status == TIMSFFI_OK && count > 0);
  }

  void BrukerTimsFile::load(const String& path, MSExperiment& exp, const Config& config)
  {
    // Placeholder — will be replaced in Task 8 with full dispatch logic
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::transform(const String& /*path*/, Interfaces::IMSDataConsumer* /*consumer*/, const Config& /*config*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::loadDDA_(void* /*handle*/, MSExperiment& /*exp*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::loadDIA_(void* /*handle*/, MSExperiment& /*exp*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::loadFrames_(void* /*handle*/, MSExperiment& /*exp*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::frameToSpectrum_(void* /*handle*/, const void* /*frame*/, MSSpectrum& /*spec*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

} // namespace OpenMS

#endif // WITH_TIMSRUST
```

- [ ] **Step 3: Register in sources.cmake (conditional)**

In `src/openms/source/FORMAT/sources.cmake`, after the `WITH_PARQUET` block (after line 121):

```cmake
if (WITH_TIMSRUST)
  list(APPEND sources_list BrukerTimsFile.cpp)
endif()
```

In `src/openms/include/OpenMS/FORMAT/sources.cmake`, after the `WITH_PARQUET` block (after line 133):

```cmake
if (WITH_TIMSRUST)
  list(APPEND sources_list_h BrukerTimsFile.h)
endif()
```

- [ ] **Step 4: Update FileHandler dispatch to use BrukerTimsFile**

In `src/openms/source/FORMAT/FileHandler.cpp`, replace the `BRUKER_TDF` stub from Task 3 with:

```cpp
#ifdef WITH_TIMSRUST
      case FileTypes::BRUKER_TDF:
      {
        BrukerTimsFile f;
        f.setLogType(log);
        f.load(filename, exp);
      }
      break;
#endif
```

Add the include at the top of FileHandler.cpp:

```cpp
#ifdef WITH_TIMSRUST
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif
```

- [ ] **Step 5: Verify build (WITH_TIMSRUST=OFF)**

Run: `cmake -B OpenMS-build -DWITH_TIMSRUST=OFF && cmake --build OpenMS-build -j$(nproc)`
Expected: Builds successfully. BrukerTimsFile is not compiled. FileHandler's BRUKER_TDF case is excluded.

- [ ] **Step 6: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h \
  src/openms/source/FORMAT/BrukerTimsFile.cpp \
  src/openms/source/FORMAT/sources.cmake \
  src/openms/include/OpenMS/FORMAT/sources.cmake \
  src/openms/source/FORMAT/FileHandler.cpp
git commit -m "feat: add BrukerTimsFile skeleton with RAII wrappers and FileHandler dispatch"
```

---

### Task 5: Frame-to-Spectrum Conversion (frameToSpectrum_)

This is the core building block used by all three load paths (DDA MS1, DIA, and raw frames).

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`
- Create: `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`
- Modify: `src/tests/class_tests/openms/executables.cmake`

- [ ] **Step 1: Create test file skeleton**

Create `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;
using namespace std;

START_TEST(BrukerTimsFile, "$Id$")

START_SECTION(void load(const String& path, MSExperiment& exp, const Config& config))
{
  // Test: invalid path throws FileNotReadable
  BrukerTimsFile f;
  MSExperiment exp;
  TEST_EXCEPTION(Exception::FileNotReadable, f.load("/nonexistent/path.d", exp));
}
END_SECTION

START_SECTION([FileHandler] BRUKER_TDF detection)
{
  // Test: .d suffix is detected as BRUKER_TDF
  TEST_EQUAL(FileHandler::getTypeByFileName("sample.d"), FileTypes::BRUKER_TDF);
  TEST_EQUAL(FileHandler::getTypeByFileName("/path/to/experiment.d"), FileTypes::BRUKER_TDF);

  // Test: non-.d suffixes are not BRUKER_TDF
  TEST_NOT_EQUAL(FileHandler::getTypeByFileName("sample.mzML"), FileTypes::BRUKER_TDF);
}
END_SECTION

END_TEST

#else // WITH_TIMSRUST

// Minimal test when timsrust is not available
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

START_TEST(BrukerTimsFile, "$Id$")
// No tests when WITH_TIMSRUST is off
END_TEST

#endif // WITH_TIMSRUST
```

- [ ] **Step 2: Register test in executables.cmake**

In `src/tests/class_tests/openms/executables.cmake`, find the `format_executables_list` section and add alphabetically:

```cmake
BrukerTimsFile_test
```

- [ ] **Step 3: Implement frameToSpectrum_**

In `src/openms/source/FORMAT/BrukerTimsFile.cpp`, replace the `frameToSpectrum_` stub:

```cpp
  void BrukerTimsFile::frameToSpectrum_(void* handle, const void* frame_ptr, MSSpectrum& spec) const
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);
    const tims_frame* frame = static_cast<const tims_frame*>(frame_ptr);

    spec.clear(true);
    spec.setRT(frame->rt_seconds);
    spec.setMSLevel(frame->ms_level);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);

    if (frame->num_peaks == 0) return;

    // Batch-convert TOF indices to m/z (check return value!)
    std::vector<double> mz_values(frame->num_peaks);
    timsffi_status mz_status = tims_convert_tof_to_mz_array(ds, frame->tof_indices, frame->num_peaks, mz_values.data());
    if (mz_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "TOF to m/z conversion failed: " + getTimsError(ds));

    // Batch-convert scan indices to IM (use batch API for efficiency)
    std::vector<uint32_t> scan_indices(frame->num_scans);
    std::iota(scan_indices.begin(), scan_indices.end(), 0u);
    std::vector<double> scan_im(frame->num_scans);
    timsffi_status im_status = tims_convert_scan_to_im_array(ds, scan_indices.data(), frame->num_scans, scan_im.data());
    if (im_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Scan to IM conversion failed: " + getTimsError(ds));

    // Build per-peak IM values from scan offsets
    std::vector<float> im_values(frame->num_peaks);
    for (uint32_t scan_idx = 0; scan_idx < frame->num_scans; ++scan_idx)
    {
      uint64_t start = frame->scan_offsets[scan_idx];
      uint64_t end = frame->scan_offsets[scan_idx + 1];
      for (uint64_t p = start; p < end; ++p)
      {
        im_values[p] = static_cast<float>(scan_im[scan_idx]);
      }
    }

    // Fill spectrum peaks
    spec.reserve(frame->num_peaks);
    for (uint32_t i = 0; i < frame->num_peaks; ++i)
    {
      Peak1D peak;
      peak.setMZ(mz_values[i]);
      peak.setIntensity(static_cast<double>(frame->intensities[i]));
      spec.push_back(peak);
    }

    // Set IM float data array with correct CV term
    DataArrays::FloatDataArray im_array;
    im_array.resize(frame->num_peaks);
    std::copy(im_values.begin(), im_values.end(), im_array.begin());
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    spec.getFloatDataArrays().push_back(std::move(im_array));
  }
```

- [ ] **Step 4: Build**

Run: `cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc) 2>&1 | tail -10`
Expected: Builds if WITH_TIMSRUST=ON and library is available, or compiles the empty test stub if OFF.

- [ ] **Step 5: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp \
  src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp \
  src/tests/class_tests/openms/executables.cmake
git commit -m "feat: implement frameToSpectrum_ core conversion with RAII and IM arrays"
```

---

### Task 6: DDA Loading (loadDDA_)

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Implement loadDDA_**

Replace the `loadDDA_` stub in `BrukerTimsFile.cpp`:

```cpp
  void BrukerTimsFile::loadDDA_(void* handle, MSExperiment& exp)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // --- MS1 frames (CONCATENATED format) ---
    unsigned int ms1_count = 0;
    tims_frame* ms1_frames = nullptr;
    timsffi_status status = tims_get_frames_by_level(ds, 1, &ms1_count, &ms1_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS1 frames: " + getTimsError(ds));

    // RAII guard for ms1_frames
    auto ms1_guard = [&]() { if (ms1_frames) tims_free_frame_array(ds, ms1_frames, ms1_count); };
    struct Ms1Guard { decltype(ms1_guard)& f; ~Ms1Guard() { f(); } } ms1_g{ms1_guard};

    exp.reserveSpaceSpectra(ms1_count + tims_num_spectra(ds));

    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      frameToSpectrum_(handle, &ms1_frames[i], spec);
      exp.addSpectrum(std::move(spec));
    }

    // --- MS2 spectra (one per precursor, scalar IM) ---
    unsigned int num_ms2 = tims_num_spectra(ds);
    for (unsigned int i = 0; i < num_ms2; ++i)
    {
      tims_spectrum ts;
      status = tims_get_spectrum(ds, i, &ts);
      if (status != TIMSFFI_OK)
      {
        OPENMS_LOG_WARN << "Warning: failed to read MS2 spectrum " << i << ": " << getTimsError(ds) << std::endl;
        continue;
      }

      MSSpectrum spec;
      spec.setRT(ts.rt_seconds);
      spec.setMSLevel(2);
      spec.setDriftTime(ts.im);
      spec.setDriftTimeUnit(DriftTimeUnit::VSSC);

      // Copy peak data (float -> double widening)
      spec.reserve(ts.num_peaks);
      for (uint32_t p = 0; p < ts.num_peaks; ++p)
      {
        Peak1D peak;
        peak.setMZ(static_cast<double>(ts.mz[p]));
        peak.setIntensity(static_cast<double>(ts.intensity[p]));
        spec.push_back(peak);
      }

      // Precursor metadata
      // Note: isolation_mz is the quadrupole center; precursor_mz is the selected ion.
      // mzML Precursor::setMZ is the isolation target, offsets are relative to it.
      Precursor prec;
      prec.setMZ(ts.isolation_mz);  // isolation window center
      prec.setCharge(static_cast<int>(ts.charge));
      if (!std::isnan(ts.precursor_intensity))
        prec.setIntensity(static_cast<float>(ts.precursor_intensity));
      prec.setDriftTime(ts.im);
      prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
      prec.setIsolationWindowLowerOffset(ts.isolation_width / 2.0);
      prec.setIsolationWindowUpperOffset(ts.isolation_width / 2.0);
      // Store selected ion m/z as user param if different from isolation center
      if (std::abs(ts.precursor_mz - ts.isolation_mz) > 1e-6)
        prec.setMetaValue("selected_ion_mz", ts.precursor_mz);

      std::vector<Precursor> precursors;
      precursors.push_back(std::move(prec));
      spec.setPrecursors(std::move(precursors));

      exp.addSpectrum(std::move(spec));
    }
  }
```

- [ ] **Step 2: Build**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: Compiles successfully.

- [ ] **Step 3: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat: implement DDA-PASEF loading (MS1 CONCATENATED + MS2 spectrum-level)"
```

---

### Task 7: DIA Loading (loadDIA_)

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Implement loadDIA_**

Replace the `loadDIA_` stub:

```cpp
  void BrukerTimsFile::loadDIA_(void* handle, MSExperiment& exp)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // --- MS1 frames (same as DDA) ---
    unsigned int ms1_count = 0;
    tims_frame* ms1_frames = nullptr;
    timsffi_status status = tims_get_frames_by_level(ds, 1, &ms1_count, &ms1_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS1 frames: " + getTimsError(ds));
    auto ms1_guard = [&]() { if (ms1_frames) tims_free_frame_array(ds, ms1_frames, ms1_count); };
    struct Ms1Guard { decltype(ms1_guard)& f; ~Ms1Guard() { f(); } } ms1_g{ms1_guard};

    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      frameToSpectrum_(handle, &ms1_frames[i], spec);
      exp.addSpectrum(std::move(spec));
    }

    // --- SWATH windows ---
    unsigned int win_count = 0;
    tims_swath_window* windows = nullptr;
    status = tims_get_swath_windows(ds, &win_count, &windows);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read SWATH windows: " + getTimsError(ds));
    auto win_guard = [&]() { if (windows) tims_free_swath_windows(ds, windows); };
    struct WinGuard { decltype(win_guard)& f; ~WinGuard() { f(); } } win_g{win_guard};

    // --- MS2 frames (split by SWATH window) ---
    unsigned int ms2_count = 0;
    tims_frame* ms2_frames = nullptr;
    status = tims_get_frames_by_level(ds, 2, &ms2_count, &ms2_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS2 frames: " + getTimsError(ds));
    auto ms2_guard = [&]() { if (ms2_frames) tims_free_frame_array(ds, ms2_frames, ms2_count); };
    struct Ms2Guard { decltype(ms2_guard)& f; ~Ms2Guard() { f(); } } ms2_g{ms2_guard};

    for (unsigned int fi = 0; fi < ms2_count; ++fi)
    {
      const tims_frame& frame = ms2_frames[fi];
      if (frame.num_peaks == 0) continue;

      // Convert TOF -> m/z for entire frame (check return value!)
      std::vector<double> mz_values(frame.num_peaks);
      timsffi_status mz_status = tims_convert_tof_to_mz_array(ds, frame.tof_indices, frame.num_peaks, mz_values.data());
      if (mz_status != TIMSFFI_OK)
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "TOF to m/z conversion failed: " + getTimsError(ds));

      // Batch scan -> IM conversion
      std::vector<uint32_t> scan_indices(frame.num_scans);
      std::iota(scan_indices.begin(), scan_indices.end(), 0u);
      std::vector<double> scan_im(frame.num_scans);
      timsffi_status im_status = tims_convert_scan_to_im_array(ds, scan_indices.data(), frame.num_scans, scan_im.data());
      if (im_status != TIMSFFI_OK)
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "Scan to IM conversion failed: " + getTimsError(ds));

      // Build per-peak IM from scan offsets
      std::vector<double> im_values(frame.num_peaks);
      for (uint32_t scan_idx = 0; scan_idx < frame.num_scans; ++scan_idx)
      {
        uint64_t start = frame.scan_offsets[scan_idx];
        uint64_t end = frame.scan_offsets[scan_idx + 1];
        for (uint64_t p = start; p < end; ++p)
        {
          im_values[p] = scan_im[scan_idx];
        }
      }

      // Split peaks by SWATH window IM bounds
      for (unsigned int wi = 0; wi < win_count; ++wi)
      {
        if (windows[wi].is_ms1) continue;

        MSSpectrum spec;
        spec.setRT(frame.rt_seconds);
        spec.setMSLevel(2);
        spec.setDriftTimeUnit(DriftTimeUnit::VSSC);

        // Set isolation window
        Precursor prec;
        prec.setMZ(windows[wi].mz_center);
        prec.setIsolationWindowLowerOffset(windows[wi].mz_center - windows[wi].mz_lower);
        prec.setIsolationWindowUpperOffset(windows[wi].mz_upper - windows[wi].mz_center);
        spec.setPrecursors({prec});

        DataArrays::FloatDataArray im_array;
        IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);

        for (uint32_t p = 0; p < frame.num_peaks; ++p)
        {
          // Check if peak falls within this window's IM bounds
          if (im_values[p] >= windows[wi].im_lower && im_values[p] <= windows[wi].im_upper)
          {
            Peak1D peak;
            peak.setMZ(mz_values[p]);
            peak.setIntensity(static_cast<double>(frame.intensities[p]));
            spec.push_back(peak);
            im_array.push_back(static_cast<float>(im_values[p]));
          }
        }

        if (!spec.empty())
        {
          spec.getFloatDataArrays().push_back(std::move(im_array));
          exp.addSpectrum(std::move(spec));
        }
      }
    }
  }
```

- [ ] **Step 2: Build**

Run: `cmake --build OpenMS-build -j$(nproc)`

- [ ] **Step 3: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat: implement DIA-PASEF loading with SWATH window splitting and per-peak IM"
```

---

### Task 8: Raw Frame Export and load() Orchestration

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Implement loadFrames_**

Replace the `loadFrames_` stub:

```cpp
  void BrukerTimsFile::loadFrames_(void* handle, MSExperiment& exp)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // Load all frames for each MS level
    for (uint8_t level = 1; level <= 2; ++level)
    {
      unsigned int count = 0;
      tims_frame* frames = nullptr;
      timsffi_status status = tims_get_frames_by_level(ds, level, &count, &frames);
      if (status != TIMSFFI_OK)
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "Failed to read frames at level " + String(level) + ": " + getTimsError(ds));
      auto guard = [&]() { if (frames) tims_free_frame_array(ds, frames, count); };
      struct Guard { decltype(guard)& f; ~Guard() { f(); } } g{guard};

      for (unsigned int i = 0; i < count; ++i)
      {
        MSSpectrum spec;
        frameToSpectrum_(handle, &frames[i], spec);
        exp.addSpectrum(std::move(spec));
      }
    }
  }
```

- [ ] **Step 2: Fix load() orchestration**

Replace the placeholder `load()` logic. Compute `isDIA_` once and use clean dispatch:

```cpp
  void BrukerTimsFile::load(const String& path, MSExperiment& exp, const Config& config)
  {
    TimsDatasetHandle ds;
    TimsConfigHandle cfg;
    openDataset(path, config, ds, cfg);

    bool is_dia = isDIA_(ds.ptr);
    Config::ExportMode mode = config.export_mode;

    if (mode == Config::FRAME)
    {
      loadFrames_(ds.ptr, exp);
    }
    else if (mode == Config::SPECTRUM || (mode == Config::AUTO && !is_dia))
    {
      loadDDA_(ds.ptr, exp);  // DDA path (also used for SPECTRUM mode on DIA data)
    }
    else // AUTO + DIA
    {
      loadDIA_(ds.ptr, exp);
    }

    // Sort by RT, interleaved across MS levels
    exp.sortSpectra(true);
  }
```

- [ ] **Step 3: Build**

Run: `cmake --build OpenMS-build -j$(nproc)`

- [ ] **Step 4: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat: implement loadFrames_ and fix load() AUTO/SPECTRUM/FRAME dispatch"
```

---

## Chunk 3: Streaming, FileConverter, and Tests

### Task 9: Streaming via transform()

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp`

- [ ] **Step 1: Implement transform()**

Replace the `transform()` stub. The implementation follows the same logic as `load()` but feeds spectra to a consumer instead of accumulating them:

```cpp
  void BrukerTimsFile::transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config)
  {
    TimsDatasetHandle ds;
    TimsConfigHandle cfg;
    openDataset(path, config, ds, cfg);

    tims_file_info_t info;
    timsffi_status status = tims_file_info(ds.ptr, &info);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read file info: " + getTimsError(ds.ptr));

    bool is_dia = isDIA_(ds.ptr);
    Config::ExportMode mode = config.export_mode;

    // Compute expected size
    size_t expected = 0;
    if (mode == Config::FRAME)
    {
      expected = info.num_frames;
    }
    else if (is_dia && mode != Config::SPECTRUM)
    {
      // DIA: MS1 frames + MS2 frames * windows
      unsigned int win_count = 0;
      tims_swath_window* windows = nullptr;
      tims_get_swath_windows(ds.ptr, &win_count, &windows);
      if (windows) tims_free_swath_windows(ds.ptr, windows);
      expected = info.ms1.count + info.ms2.count * win_count;
    }
    else
    {
      expected = info.ms1.count + info.num_spectra_ms2;
    }

    consumer->setExpectedSize(expected, 0);
    consumer->setExperimentalSettings(ExperimentalSettings());

    // NOTE: This loads into a temporary experiment then feeds to consumer.
    // Not truly constant-memory — a future optimization should iterate
    // frame-by-frame and call consumer->consumeSpectrum() inline.
    OPENMS_LOG_INFO << "BrukerTimsFile::transform(): loading full dataset (streaming optimization pending)" << std::endl;
    MSExperiment exp;
    if (mode == Config::FRAME)
      loadFrames_(ds.ptr, exp);
    else if (is_dia && mode != Config::SPECTRUM)
      loadDIA_(ds.ptr, exp);
    else
      loadDDA_(ds.ptr, exp);

    exp.sortSpectra(true);
    for (auto& spec : exp)
    {
      consumer->consumeSpectrum(spec);
    }
  }
```

- [ ] **Step 2: Build**

Run: `cmake --build OpenMS-build -j$(nproc)`

- [ ] **Step 3: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat: implement transform() streaming via IMSDataConsumer"
```

---

### Task 10: FileConverter Integration

**Files:**
- Modify: `src/topp/FileConverter.cpp:152` (input formats)
- Modify: `src/topp/FileConverter.cpp:148-180` (parameter registration)
- Modify: `src/topp/FileConverter.cpp:394-456` (low-memory branch)

- [ ] **Step 1: Add "d" to input formats**

In `src/topp/FileConverter.cpp`, line 152, add `"d"` to the `input_formats` vector:

```cpp
    vector<String> input_formats = {"mzML", "mzXML", "mgf", "msp", "raw", "cachedMzML", "mzData", "dta", "dta2d", "featureXML", "consensusXML", "ms2", "fid", "d", "tsv", "peplist", "kroenik", "edta", "oms", "sqMass"};
```

- [ ] **Step 2: Add timsrust parameters (guarded)**

After the existing parameter registration block (after line ~176, before `registerTOPPSubsection_("RawToMzML", ...)`), add:

```cpp
#ifdef WITH_TIMSRUST
    registerTOPPSubsection_("timsrust", "Options for reading Bruker TimsTOF .d files (requires WITH_TIMSRUST)");
    registerIntOption_("timsrust:smoothing_window", "<int>", 0, "Smoothing window bin count (0 = library default)", false, true);
    setMinInt_("timsrust:smoothing_window", 0);
    registerIntOption_("timsrust:centroiding_window", "<int>", 0, "Centroiding window bin count (0 = library default)", false, true);
    setMinInt_("timsrust:centroiding_window", 0);
    registerDoubleOption_("timsrust:calibration_tolerance", "<float>", 0.0, "m/z calibration tolerance (0 = library default)", false, true);
    setMinFloat_("timsrust:calibration_tolerance", 0.0);
    registerStringOption_("timsrust:calibrate", "<toggle>", "false", "Enable m/z recalibration (may fail on some datasets)", false, true);
    setValidStrings_("timsrust:calibrate", {"true", "false"});
    registerStringOption_("timsrust:export_mode", "<mode>", "auto", "Export mode: auto (detect DDA/DIA), spectrum (force spectrum-level), frame (raw 4D frames)", false, true);
    setValidStrings_("timsrust:export_mode", {"auto", "spectrum", "frame"});
#endif
```

- [ ] **Step 3: Add includes**

At the top of `FileConverter.cpp`, add:

```cpp
#ifdef WITH_TIMSRUST
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif
```

- [ ] **Step 4: Extend low-memory branch for BRUKER_TDF**

In the `process_lowmemory` block (around line 401), add a new branch before the `else` (line 451):

```cpp
#ifdef WITH_TIMSRUST
      else if (in_type == FileTypes::BRUKER_TDF && out_type == FileTypes::MZML)
      {
        PlainMSDataWritingConsumer consumer(out);
        consumer.getOptions().setWriteIndex(write_scan_index);
        if (lossy_compression)
        {
          consumer.getOptions().setNumpressConfigurationMassTime(npconfig_mz);
          consumer.getOptions().setNumpressConfigurationIntensity(npconfig_int);
          consumer.getOptions().setNumpressConfigurationFloatDataArray(npconfig_fda);
          consumer.getOptions().setCompression(true);
        }
        consumer.addDataProcessing(getProcessingInfo_(DataProcessing::CONVERSION_MZML));

        BrukerTimsFile tims_file;
        tims_file.setLogType(log_type_);
        BrukerTimsFile::Config tims_config;
        tims_config.smoothing_window = static_cast<uint32_t>(getIntOption_("timsrust:smoothing_window"));
        tims_config.centroiding_window = static_cast<uint32_t>(getIntOption_("timsrust:centroiding_window"));
        tims_config.calibration_tolerance = getDoubleOption_("timsrust:calibration_tolerance");
        tims_config.calibrate = (getStringOption_("timsrust:calibrate") == "true");
        String mode = getStringOption_("timsrust:export_mode");
        if (mode == "spectrum") tims_config.export_mode = BrukerTimsFile::Config::SPECTRUM;
        else if (mode == "frame") tims_config.export_mode = BrukerTimsFile::Config::FRAME;
        else tims_config.export_mode = BrukerTimsFile::Config::AUTO;

        tims_file.transform(in, &consumer, tims_config);
        return EXECUTION_OK;
      }
#endif
```

- [ ] **Step 5: Add timsrust config to normal processing branch**

In the normal processing branch (around line 457-460), before `fh.loadExperiment(...)`, add config extraction for `.d` input:

```cpp
#ifdef WITH_TIMSRUST
    else if (in_type == FileTypes::BRUKER_TDF)
    {
      BrukerTimsFile tims_file;
      tims_file.setLogType(log_type_);
      BrukerTimsFile::Config tims_config;
      tims_config.smoothing_window = static_cast<uint32_t>(getIntOption_("timsrust:smoothing_window"));
      tims_config.centroiding_window = static_cast<uint32_t>(getIntOption_("timsrust:centroiding_window"));
      tims_config.calibration_tolerance = getDoubleOption_("timsrust:calibration_tolerance");
      tims_config.calibrate = (getStringOption_("timsrust:calibrate") == "true");
      String mode = getStringOption_("timsrust:export_mode");
      if (mode == "spectrum") tims_config.export_mode = BrukerTimsFile::Config::SPECTRUM;
      else if (mode == "frame") tims_config.export_mode = BrukerTimsFile::Config::FRAME;
      else tims_config.export_mode = BrukerTimsFile::Config::AUTO;

      tims_file.load(in, exp, tims_config);
    }
#endif
```

This goes between the `process_lowmemory` block and the generic `fh.loadExperiment()` call.

- [ ] **Step 6: Build**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: Compiles (WITH_TIMSRUST=OFF: guarded blocks excluded; WITH_TIMSRUST=ON: needs library).

- [ ] **Step 7: Commit**

```bash
git add src/topp/FileConverter.cpp
git commit -m "feat: integrate BrukerTimsFile into FileConverter with timsrust parameters"
```

---

### Task 11: Integration Test Infrastructure

**Files:**
- Modify: `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`
- Modify: `CMakeLists.txt` (test data FetchContent, after `ENABLE_TIMSRUST_TESTS`)

- [ ] **Step 1: Add test data download in CMake**

In the top-level `CMakeLists.txt`, after the `ENABLE_TIMSRUST_TESTS` option, add:

```cmake
if (ENABLE_TIMSRUST_TESTS AND WITH_TIMSRUST)
  include(FetchContent)
  set(TIMSRUST_TESTDATA_URL "https://github.com/OpenMS/timsrust_cpp_bridge/releases/download/test-data-v1")
  FetchContent_Declare(
    timsrust_testdata_dda
    URL "${TIMSRUST_TESTDATA_URL}/DDA_HeLa_50ng_5_6min.d.zip"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_Declare(
    timsrust_testdata_dia
    URL "${TIMSRUST_TESTDATA_URL}/DIA_HeLa_50ng_5_6min.d.zip"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(timsrust_testdata_dda timsrust_testdata_dia)
  set(TIMSRUST_DDA_TEST_DATA "${timsrust_testdata_dda_SOURCE_DIR}" CACHE PATH "Path to DDA test .d directory")
  set(TIMSRUST_DIA_TEST_DATA "${timsrust_testdata_dia_SOURCE_DIR}" CACHE PATH "Path to DIA test .d directory")
endif()
```

Also, in `src/tests/class_tests/openms/executables.cmake` or the test target setup, pass the paths as compile definitions so the test can access them:

```cmake
if (TARGET BrukerTimsFile_test AND TIMSRUST_DDA_TEST_DATA)
  target_compile_definitions(BrukerTimsFile_test PRIVATE
    TIMSRUST_DDA_TEST_DATA="${TIMSRUST_DDA_TEST_DATA}"
    TIMSRUST_DIA_TEST_DATA="${TIMSRUST_DIA_TEST_DATA}")
endif()
```

- [ ] **Step 2: Add integration test sections**

In `BrukerTimsFile_test.cpp`, add integration test sections guarded by a macro that checks whether test data is available. The test data paths would be passed via CMake defines. For now, add placeholder sections:

```cpp
// Integration tests (only run when ENABLE_TIMSRUST_TESTS is ON and data is available)
#ifdef TIMSRUST_DDA_TEST_DATA

START_SECTION(DDA loading integration test)
{
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(TIMSRUST_DDA_TEST_DATA, exp);

  // Verify we got spectra
  TEST_NOT_EQUAL(exp.size(), 0);

  // Check MS levels present
  bool has_ms1 = false, has_ms2 = false;
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 1) has_ms1 = true;
    if (spec.getMSLevel() == 2) has_ms2 = true;
  }
  TEST_EQUAL(has_ms1, true);
  TEST_EQUAL(has_ms2, true);

  // Check MS1 spectra have IM data in CONCATENATED format
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 1 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      break;
    }
  }

  // Check MS2 spectra have precursor and drift time
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty())
    {
      TEST_NOT_EQUAL(spec.getDriftTime(), 0.0);
      TEST_EQUAL(spec.getDriftTimeUnit(), DriftTimeUnit::VSSC);
      TEST_NOT_EQUAL(spec.getPrecursors()[0].getMZ(), 0.0);
      break;
    }
  }
}
END_SECTION

#endif // TIMSRUST_DDA_TEST_DATA

#ifdef TIMSRUST_DIA_TEST_DATA

START_SECTION(DIA loading integration test)
{
  BrukerTimsFile f;
  MSExperiment exp;
  f.load(TIMSRUST_DIA_TEST_DATA, exp);

  TEST_NOT_EQUAL(exp.size(), 0);

  // Check MS2 spectra have per-peak IM (CONCATENATED) and isolation windows
  for (const auto& spec : exp)
  {
    if (spec.getMSLevel() == 2 && !spec.empty())
    {
      TEST_EQUAL(spec.containsIMData(), true);
      TEST_NOT_EQUAL(spec.getPrecursors().size(), 0);
      break;
    }
  }
}
END_SECTION

#endif // TIMSRUST_DIA_TEST_DATA
```

- [ ] **Step 3: Commit**

```bash
git add CMakeLists.txt \
  src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp
git commit -m "test: add integration test infrastructure for Bruker TimsTOF with real data"
```

---

## Summary

| Task | Description | Key Deliverable |
|------|-------------|-----------------|
| 1 | CMake option + FetchContent + config.h.in | Build infrastructure |
| 2 | BRUKER_TDF file type registration | FileTypes enum + bindings + test |
| 3 | FileHandler directory detection + dispatch | getType() + loadExperiment() |
| 4 | BrukerTimsFile skeleton + RAII wrappers | Header, source, sources.cmake |
| 5 | frameToSpectrum_ core conversion | TOF→m/z, scan→IM, FloatDataArray |
| 6 | loadDDA_ implementation | MS1 CONCATENATED + MS2 scalar IM |
| 7 | loadDIA_ implementation | SWATH window splitting + per-peak IM |
| 8 | loadFrames_ + load() orchestration | Frame export + AUTO detection |
| 9 | transform() streaming | IMSDataConsumer integration |
| 10 | FileConverter integration | Parameters + low-memory branch |
| 11 | Integration test infrastructure | Test data download + DDA/DIA tests |
