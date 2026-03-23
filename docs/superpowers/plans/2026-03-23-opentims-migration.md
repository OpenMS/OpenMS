# opentims Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the timsrust Rust bridge with opentims (C++) plus open-source calibration converters for reading Bruker TimsTOF `.d` files.

**Architecture:** opentims handles TDF binary I/O (ZSTD decompression, memory-mapped frame reading). Two new C++ converter classes implement timsrust's linear-in-sqrt calibration math, registered via opentims's factory pattern. BrukerTimsFile.cpp is rewritten against opentims's `TimsDataHandle`/`TimsFrame` API, with DDA precursor grouping and DIA window parsing done via direct SQL queries against `analysis.tdf`.

**Tech Stack:** C++17, CMake FetchContent, opentims (C++20), SQLiteCpp (existing in OpenMS), ZSTD (bundled in opentims)

**Spec:** `docs/superpowers/specs/2026-03-23-opentims-migration-design.md`

---

## File Structure

| File | Responsibility |
|---|---|
| `cmake/cmake_findExternalLibs.cmake` | FetchContent opentims, build opentims sources into OpenMS |
| `CMakeLists.txt` (top-level) | `WITH_OPENTIMS` option, test data fetch |
| `cmake/OpenMSConfig.cmake.in` | Export `WITH_OPENTIMS` to downstream |
| `src/openms/CMakeLists.txt` | Link opentims, compile definitions |
| `src/openms/include/OpenMS/FORMAT/sources.cmake` | Register new headers |
| `src/openms/source/FORMAT/sources.cmake` | Register new sources |
| `src/openms/include/OpenMS/FORMAT/OpenTimsCalibration.h` | **NEW**: converter + factory declarations |
| `src/openms/source/FORMAT/OpenTimsCalibration.cpp` | **NEW**: calibration math implementation |
| `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h` | Updated header (remove timsrust refs, update Config) |
| `src/openms/source/FORMAT/BrukerTimsFile.cpp` | Full rewrite against opentims API |
| `src/openms/source/FORMAT/FileHandler.cpp` | `#ifdef` rename |
| `src/topp/FileConverter.cpp` | TOPP parameter rename `timsrust:` → `bruker:` |
| `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp` | Updated tests + new converter/factory tests |
| `src/tests/class_tests/openms/CMakeLists.txt` | Rename `TIMSRUST_*` compile definitions to `OPENTIMS_*` |
| `src/tests/topp/CMakeLists.txt` | Variable renames |

---

## Task 1: Generate Golden-File References

Before changing anything, capture the current timsrust-based output for regression testing.

**Files:**
- Read: `src/openms/source/FORMAT/BrukerTimsFile.cpp` (existing)
- Create: `src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dda.mzML` (reference)
- Create: `src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dia.mzML` (reference)

**Prerequisite:** `ENABLE_TIMSRUST_TESTS=ON` and test data available.

- [ ] **Step 1: Write a script to dump golden files**

Create a small C++ test or use FileConverter to load DDA and DIA test `.d` files and write them to mzML:

```bash
cd /home/sachsenb/Development/tmp/OpenMS
# Build current code
cmake --build OpenMS-build -j$(nproc)
# Generate golden files
OpenMS-build/bin/FileConverter -in "${TIMSRUST_DDA_TEST_DATA}" -out src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dda.mzML
OpenMS-build/bin/FileConverter -in "${TIMSRUST_DIA_TEST_DATA}" -out src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dia.mzML
```

Note: The `TIMSRUST_DDA_TEST_DATA` / `TIMSRUST_DIA_TEST_DATA` paths are set by CMake. Check `OpenMS-build/CMakeCache.txt` for the actual paths.

- [ ] **Step 2: Verify golden files are non-empty and contain expected data**

```bash
grep -c '<spectrum ' src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dda.mzML
grep -c '<spectrum ' src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dia.mzML
```

Both should report > 0 spectra.

- [ ] **Step 3: Commit golden files**

```bash
git add src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dda.mzML \
        src/tests/class_tests/openms/test_data/BrukerTimsFile_golden_dia.mzML
git commit -m "test: add golden-file references for opentims migration regression testing"
```

---

## Task 2: CMake — Replace timsrust with opentims FetchContent

Remove all timsrust CMake infrastructure and add opentims.

**Files:**
- Modify: `CMakeLists.txt:61-98` (options + test data)
- Modify: `cmake/cmake_findExternalLibs.cmake:276-319` (FetchContent)
- Modify: `cmake/OpenMSConfig.cmake.in:72`
- Modify: `src/openms/CMakeLists.txt:101-103,175-177` (link + compile def)
- Modify: `src/openms/include/OpenMS/FORMAT/sources.cmake:137-139`
- Modify: `src/openms/source/FORMAT/sources.cmake:126-128`

- [ ] **Step 1: Update top-level CMakeLists.txt options**

In `CMakeLists.txt`, replace lines 61-62:

```cmake
# OLD:
option(WITH_TIMSRUST "Build with Bruker TimsTOF .d file support via timsrust" ON)
option(ENABLE_TIMSRUST_TESTS "Download Bruker test data for timsrust integration tests" OFF)

# NEW:
option(WITH_OPENTIMS "Build with Bruker TimsTOF .d file support via opentims" ON)
option(ENABLE_OPENTIMS_TESTS "Download Bruker test data for opentims integration tests" OFF)
```

Replace lines 76-98 (test data fetch block):

```cmake
# OLD: if (ENABLE_TIMSRUST_TESTS AND WITH_TIMSRUST) ... endif()
# NEW:
if (ENABLE_OPENTIMS_TESTS AND WITH_OPENTIMS)
  include(FetchContent)
  set(OPENTIMS_TESTDATA_URL "https://github.com/OpenMS/timsrust_cpp_bridge/releases/download/test-data-v1")
  FetchContent_Declare(
    opentims_testdata_dda
    URL "${OPENTIMS_TESTDATA_URL}/DDA_HeLa_50ng_5_6min.d.zip"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_Declare(
    opentims_testdata_dia
    URL "${OPENTIMS_TESTDATA_URL}/DIA_HeLa_50ng_5_6min.d.zip"
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_Declare(
    opentims_testdata_fasta
    URL "${OPENTIMS_TESTDATA_URL}/human_sp.fasta"
    DOWNLOAD_NO_EXTRACT TRUE
  )
  FetchContent_MakeAvailable(opentims_testdata_dda opentims_testdata_dia opentims_testdata_fasta)
  set(OPENTIMS_DDA_TEST_DATA "${opentims_testdata_dda_SOURCE_DIR}" CACHE PATH "Path to DDA test .d directory")
  set(OPENTIMS_DIA_TEST_DATA "${opentims_testdata_dia_SOURCE_DIR}" CACHE PATH "Path to DIA test .d directory")
  set(OPENTIMS_TEST_FASTA "${opentims_testdata_fasta_SOURCE_DIR}/human_sp.fasta" CACHE PATH "Path to human Swiss-Prot FASTA for search tests")
endif()
```

- [ ] **Step 2: Replace timsrust FetchContent with opentims in cmake_findExternalLibs.cmake**

Replace lines 276-319 entirely:

```cmake
#------------------------------------------------------------------------------
# opentims (Bruker TimsTOF .d file reading)
if (WITH_OPENTIMS)
  include(FetchContent)
  # Pin to a specific commit for reproducible builds.
  # TODO: Update this hash after license clarification with opentims maintainer.
  set(OPENTIMS_COMMIT "PLACEHOLDER_COMMIT_HASH" CACHE STRING "opentims commit hash")

  FetchContent_Declare(
    opentims
    GIT_REPOSITORY https://github.com/michalsta/opentims.git
    GIT_TAG ${OPENTIMS_COMMIT}
  )
  FetchContent_GetProperties(opentims)
  if(NOT opentims_POPULATED)
    FetchContent_Populate(opentims)
  endif()

  # Build opentims C++ sources as a static library, excluding bundled sqlite
  # (we use OpenMS's existing SQLiteCpp/sqlite3 instead).
  set(_OPENTIMS_SRC_DIR "${opentims_SOURCE_DIR}/src/opentims++")
  # Hand-maintained file list (not GLOB) per spec. Update when upgrading opentims.
  # Verify the actual filenames against the pinned commit at implementation time.
  set(_OPENTIMS_SOURCES
    "${_OPENTIMS_SRC_DIR}/opentims.cpp"
    "${_OPENTIMS_SRC_DIR}/tof2mz_converter.cpp"
    "${_OPENTIMS_SRC_DIR}/scan2inv_ion_mobility_converter.cpp"
    "${_OPENTIMS_SRC_DIR}/converters.cpp"
    "${_OPENTIMS_SRC_DIR}/so_manager.cpp"
    # Add other .cpp files from opentims as needed; exclude sqlite3.c and pybind11
  )

  add_library(opentims_cpp STATIC ${_OPENTIMS_SOURCES})
  target_include_directories(opentims_cpp PUBLIC "${_OPENTIMS_SRC_DIR}")
  # Use OpenMS's sqlite3 headers instead of opentims's bundled copy
  target_include_directories(opentims_cpp PRIVATE
    "${CMAKE_SOURCE_DIR}/src/openms/extern/SQLiteCpp/sqlite3")
  # opentims uses C++20 features. PRIVATE to avoid forcing C++20 on all of OpenMS.
  # If opentims headers used by OpenMS code require C++20, change to PUBLIC and
  # verify the full OpenMS build still works.
  target_compile_features(opentims_cpp PRIVATE cxx_std_20)
  # Link against OpenMS's sqlite3 static library
  target_link_libraries(opentims_cpp PRIVATE sqlite3)
  # Suppress warnings in third-party code
  target_compile_options(opentims_cpp PRIVATE -w)

  message(STATUS "Found opentims: ${opentims_SOURCE_DIR}")
endif()
```

- [ ] **Step 3: Update OpenMSConfig.cmake.in**

In `cmake/OpenMSConfig.cmake.in`, replace line 72:

```cmake
# OLD: set(WITH_TIMSRUST @WITH_TIMSRUST@)
# NEW:
set(WITH_OPENTIMS @WITH_OPENTIMS@)
```

- [ ] **Step 4: Update src/openms/CMakeLists.txt linkage**

Replace lines 101-103:

```cmake
# OLD:
if (WITH_TIMSRUST)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES timsrust_cpp_bridge::timsrust_cpp_bridge)
endif()

# NEW:
if (WITH_OPENTIMS)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES opentims_cpp)
endif()
```

Replace lines 175-177:

```cmake
# OLD:
if (WITH_TIMSRUST)
    target_compile_definitions(OpenMS PUBLIC WITH_TIMSRUST=1)
endif()

# NEW:
if (WITH_OPENTIMS)
    target_compile_definitions(OpenMS PUBLIC WITH_OPENTIMS=1)
endif()
```

- [ ] **Step 5: Update sources.cmake files**

In `src/openms/include/OpenMS/FORMAT/sources.cmake`, replace lines 137-139:

```cmake
# OLD:
if (WITH_TIMSRUST)
  list(APPEND sources_list_h BrukerTimsFile.h)
endif()

# NEW:
if (WITH_OPENTIMS)
  list(APPEND sources_list_h BrukerTimsFile.h)
  list(APPEND sources_list_h OpenTimsCalibration.h)
endif()
```

In `src/openms/source/FORMAT/sources.cmake`, replace lines 126-128:

```cmake
# OLD:
if (WITH_TIMSRUST)
  list(APPEND sources_list BrukerTimsFile.cpp)
endif()

# NEW:
if (WITH_OPENTIMS)
  list(APPEND sources_list BrukerTimsFile.cpp)
  list(APPEND sources_list OpenTimsCalibration.cpp)
endif()
```

- [ ] **Step 6: Verify CMake configures**

```bash
cd /home/sachsenb/Development/tmp/OpenMS
cmake -S . -B OpenMS-build -DWITH_OPENTIMS=OFF
```

Expected: configures successfully with opentims disabled (no fetch).

- [ ] **Step 7: Commit**

```bash
git add CMakeLists.txt cmake/cmake_findExternalLibs.cmake cmake/OpenMSConfig.cmake.in \
        src/openms/CMakeLists.txt \
        src/openms/include/OpenMS/FORMAT/sources.cmake \
        src/openms/source/FORMAT/sources.cmake
git commit -m "build: replace timsrust with opentims FetchContent infrastructure"
```

---

## Task 3: Implement Open-Source Calibration Converters

**Files:**
- Create: `src/openms/include/OpenMS/FORMAT/OpenTimsCalibration.h`
- Create: `src/openms/source/FORMAT/OpenTimsCalibration.cpp`

These classes implement opentims's `Tof2MzConverter` and `Scan2InvIonMobilityConverter` interfaces using timsrust's linear-in-sqrt calibration math. They read parameters from the `analysis.tdf` SQLite metadata.

- [ ] **Step 1: Create OpenTimsCalibration.h**

Create `src/openms/include/OpenMS/FORMAT/OpenTimsCalibration.h`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <tof2mz_converter.h>
#include <scan2inv_ion_mobility_converter.h>
#include <string>

namespace OpenMS
{

  /// Open-source TOF-to-m/z converter using linear-in-sqrt calibration (matching timsrust).
  /// Reads MzAcqRangeLower, MzAcqRangeUpper, DigitizerNumSamples from GlobalMetadata.
  class OpenSourceTof2MzConverter : public Tof2MzConverter
  {
  public:
    OpenSourceTof2MzConverter(double mz_min, double mz_max, uint32_t tof_max_index,
                              bool is_otof_control);

    void convert(uint32_t frame_id, double* mzs, const double* tofs, uint32_t size) override;
    void convert(uint32_t frame_id, double* mzs, const uint32_t* tofs, uint32_t size) override;
    void inverse_convert(uint32_t frame_id, double* tofs, const double* mzs, uint32_t size) override;
    std::string description() const override;

    /// Update calibration parameters (used by OLS recalibration)
    void updateCalibration(double new_intercept, double new_slope);
    double intercept() const { return intercept_; }
    double slope() const { return slope_; }

  private:
    double intercept_;
    double slope_;
  };

  /// Open-source scan-to-inverse-ion-mobility converter (linear model).
  class OpenSourceScan2ImConverter : public Scan2InvIonMobilityConverter
  {
  public:
    OpenSourceScan2ImConverter(double im_min, double im_max, uint32_t scan_max_index);

    void convert(uint32_t frame_id, double* inv_ion_mobilities, const double* scans, uint32_t size) override;
    void convert(uint32_t frame_id, double* inv_ion_mobilities, const uint32_t* scans, uint32_t size) override;
    void inverse_convert(uint32_t frame_id, double* scans, const double* inv_ion_mobilities, uint32_t size) override;
    std::string description() const override;

  private:
    double intercept_;
    double slope_;
  };

  /// Factory for OpenSourceTof2MzConverter. Reads calibration from analysis.tdf.
  class OpenSourceTof2MzConverterFactory : public Tof2MzConverterFactory
  {
  public:
    std::unique_ptr<Tof2MzConverter> produce(TimsDataHandle& TDH,
      pressure_compensation_strategy pcs = NoPressureCompensation) override;
  };

  /// Factory for OpenSourceScan2ImConverter. Reads calibration from analysis.tdf.
  class OpenSourceScan2ImConverterFactory : public Scan2InvIonMobilityConverterFactory
  {
  public:
    std::unique_ptr<Scan2InvIonMobilityConverter> produce(TimsDataHandle& TDH,
      pressure_compensation_strategy pcs = NoPressureCompensation) override;
  };

} // namespace OpenMS

#endif // WITH_OPENTIMS
```

Note: The exact opentims base class names and `produce()` signatures must be verified against the pinned opentims commit. Adjust if the API differs.

- [ ] **Step 2: Create OpenTimsCalibration.cpp**

Create `src/openms/source/FORMAT/OpenTimsCalibration.cpp`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/OpenTimsCalibration.h>
#include <opentims.h>
#include <SQLiteCpp/SQLiteCpp.h>
#include <cmath>
#include <stdexcept>

namespace OpenMS
{

  // --- OpenSourceTof2MzConverter ---

  OpenSourceTof2MzConverter::OpenSourceTof2MzConverter(
    double mz_min, double mz_max, uint32_t tof_max_index, bool is_otof_control)
  {
    if (is_otof_control)
    {
      mz_min -= 5.0;
      mz_max += 5.0;
    }
    intercept_ = std::sqrt(mz_min);
    slope_ = (std::sqrt(mz_max) - std::sqrt(mz_min)) / static_cast<double>(tof_max_index);
  }

  void OpenSourceTof2MzConverter::convert(uint32_t /*frame_id*/, double* mzs,
    const double* tofs, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      double val = intercept_ + slope_ * tofs[i];
      mzs[i] = val * val;
    }
  }

  void OpenSourceTof2MzConverter::convert(uint32_t /*frame_id*/, double* mzs,
    const uint32_t* tofs, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      double val = intercept_ + slope_ * static_cast<double>(tofs[i]);
      mzs[i] = val * val;
    }
  }

  void OpenSourceTof2MzConverter::inverse_convert(uint32_t /*frame_id*/, double* tofs,
    const double* mzs, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      tofs[i] = (std::sqrt(mzs[i]) - intercept_) / slope_;
    }
  }

  std::string OpenSourceTof2MzConverter::description() const
  {
    return "OpenSourceTof2MzConverter (linear-in-sqrt, OpenMS)";
  }

  void OpenSourceTof2MzConverter::updateCalibration(double new_intercept, double new_slope)
  {
    intercept_ = new_intercept;
    slope_ = new_slope;
  }

  // --- OpenSourceScan2ImConverter ---

  OpenSourceScan2ImConverter::OpenSourceScan2ImConverter(
    double im_min, double im_max, uint32_t scan_max_index)
  {
    intercept_ = im_max;
    slope_ = (im_min - im_max) / static_cast<double>(scan_max_index);
  }

  void OpenSourceScan2ImConverter::convert(uint32_t /*frame_id*/, double* inv_ion_mobilities,
    const double* scans, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = intercept_ + slope_ * scans[i];
    }
  }

  void OpenSourceScan2ImConverter::convert(uint32_t /*frame_id*/, double* inv_ion_mobilities,
    const uint32_t* scans, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = intercept_ + slope_ * static_cast<double>(scans[i]);
    }
  }

  void OpenSourceScan2ImConverter::inverse_convert(uint32_t /*frame_id*/, double* scans,
    const double* inv_ion_mobilities, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      scans[i] = (inv_ion_mobilities[i] - intercept_) / slope_;
    }
  }

  std::string OpenSourceScan2ImConverter::description() const
  {
    return "OpenSourceScan2ImConverter (linear, OpenMS)";
  }

  // --- Factories ---

  std::unique_ptr<Tof2MzConverter> OpenSourceTof2MzConverterFactory::produce(
    TimsDataHandle& TDH, pressure_compensation_strategy /*pcs*/)
  {
    std::string tdf_path = TDH.tims_dir_path + "/analysis.tdf";
    SQLite::Database db(tdf_path, SQLite::OPEN_READONLY);

    // Read calibration parameters from GlobalMetadata
    double mz_min = 0, mz_max = 0;
    uint32_t tof_max = 0;
    bool is_otof = false;

    SQLite::Statement query(db, "SELECT Key, Value FROM GlobalMetadata "
      "WHERE Key IN ('MzAcqRangeLower','MzAcqRangeUpper','DigitizerNumSamples','AcquisitionSoftware')");

    while (query.executeStep())
    {
      std::string key = query.getColumn(0).getString();
      std::string val = query.getColumn(1).getString();
      if (key == "MzAcqRangeLower") mz_min = std::stod(val);
      else if (key == "MzAcqRangeUpper") mz_max = std::stod(val);
      else if (key == "DigitizerNumSamples") tof_max = static_cast<uint32_t>(std::stoul(val));
      else if (key == "AcquisitionSoftware") is_otof = (val == "Bruker otofControl");
    }

    if (mz_min <= 0 || mz_max <= 0 || tof_max == 0)
      throw std::runtime_error("OpenSourceTof2MzConverterFactory: missing calibration metadata in " + tdf_path);

    return std::make_unique<OpenSourceTof2MzConverter>(mz_min, mz_max, tof_max, is_otof);
  }

  std::unique_ptr<Scan2InvIonMobilityConverter> OpenSourceScan2ImConverterFactory::produce(
    TimsDataHandle& TDH, pressure_compensation_strategy /*pcs*/)
  {
    std::string tdf_path = TDH.tims_dir_path + "/analysis.tdf";
    SQLite::Database db(tdf_path, SQLite::OPEN_READONLY);

    double im_min = 0, im_max = 0;

    SQLite::Statement meta_query(db, "SELECT Key, Value FROM GlobalMetadata "
      "WHERE Key IN ('OneOverK0AcqRangeLower','OneOverK0AcqRangeUpper')");
    while (meta_query.executeStep())
    {
      std::string key = meta_query.getColumn(0).getString();
      std::string val = meta_query.getColumn(1).getString();
      if (key == "OneOverK0AcqRangeLower") im_min = std::stod(val);
      else if (key == "OneOverK0AcqRangeUpper") im_max = std::stod(val);
    }

    SQLite::Statement scan_query(db, "SELECT MAX(NumScans) FROM Frames");
    uint32_t scan_max = 0;
    if (scan_query.executeStep())
    {
      scan_max = static_cast<uint32_t>(scan_query.getColumn(0).getInt());
    }

    if (im_min <= 0 || im_max <= 0 || scan_max == 0)
      throw std::runtime_error("OpenSourceScan2ImConverterFactory: missing calibration metadata in " + tdf_path);

    return std::make_unique<OpenSourceScan2ImConverter>(im_min, im_max, scan_max);
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
```

- [ ] **Step 3: Verify files compile (with opentims disabled — just syntax check)**

```bash
# This won't fully compile until opentims is actually fetched, but confirms
# no obvious syntax errors in the #ifdef-guarded code when disabled.
cmake --build OpenMS-build -j$(nproc) 2>&1 | head -50
```

- [ ] **Step 4: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/OpenTimsCalibration.h \
        src/openms/source/FORMAT/OpenTimsCalibration.cpp
git commit -m "feat: add open-source TOF-to-m/z and scan-to-IM calibration converters for opentims"
```

---

## Task 4: Update BrukerTimsFile.h

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h`

- [ ] **Step 1: Update the header**

Replace the entire file content. Key changes:
- `#ifdef WITH_TIMSRUST` → `#ifdef WITH_OPENTIMS`
- Remove `smoothing_window` and `centroiding_window` from `Config`
- Update doxygen to remove timsrust references
- Change private methods to use opentims types instead of `void*`

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <cstdint>

// Forward declaration — opentims types
class TimsDataHandle;

namespace OpenMS
{

  /**
   * @brief Reader for Bruker TimsTOF .d directories via opentims.
   *
   * Supports DDA-PASEF, DIA-PASEF, and raw frame-level 4D access.
   * Ion mobility data is stored in VSSC (1/K0) units using CONCATENATED format
   * for MS1 and DIA MS2, and scalar drift times for DDA MS2.
   *
   * @section export_modes Export modes
   *
   * - **AUTO** (default): detects DDA vs DIA by querying for SWATH windows.
   *   DDA uses per-precursor MS2 spectra with scalar drift times; DIA uses
   *   per-frame MS2 spectra split by SWATH isolation window with per-peak IM.
   * - **SPECTRUM**: forces the DDA (per-precursor) path regardless of
   *   acquisition type.
   * - **FRAME**: returns raw 4D frames as CONCATENATED spectra (per-peak IM
   *   in FloatDataArray) for both MS1 and MS2, without precursor grouping or
   *   SWATH splitting.
   */
  class OPENMS_DLLAPI BrukerTimsFile : public ProgressLogger
  {
  public:
    /// Processing and export configuration
    struct Config
    {
      double calibration_tolerance = 0.0;  ///< m/z recalibration tolerance (0 = default 0.1 Da)
      bool calibrate = false;              ///< Enable OLS m/z recalibration (off by default)

      float ms1_centroid_mz_ppm = 0.0f; ///< MS1 IM-centroiding m/z tolerance in ppm (0 = disabled, suggested: 5.0)
      float ms1_centroid_im_pct = 0.0f;  ///< MS1 IM-centroiding ion mobility tolerance in percent (0 = disabled, suggested: 3.0)

      enum ExportMode { AUTO, SPECTRUM, FRAME };
      ExportMode export_mode = AUTO;
    };

    /// Load entire .d directory into MSExperiment
    void load(const String& path, MSExperiment& exp);
    /// @overload with explicit configuration
    void load(const String& path, MSExperiment& exp, const Config& config);

    /// Feed spectra from a .d directory to a consumer.
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer);
    /// @overload with explicit configuration
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config);

  private:
    void loadDDA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);
    void loadDIA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);
    void loadFrames_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);
    bool isDIA_(const String& tdf_path) const;
  };

} // namespace OpenMS

#endif // WITH_OPENTIMS
```

- [ ] **Step 2: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h
git commit -m "refactor: update BrukerTimsFile header for opentims (remove timsrust types)"
```

---

## Task 5: Rewrite BrukerTimsFile.cpp

This is the largest task. The full rewrite replaces all timsrust C FFI calls with opentims `TimsDataHandle`/`TimsFrame` API and direct SQL queries.

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp` (full rewrite)

**Key reference:** See spec Section 4 "DDA MS2 Spectrum Reconstruction" for the detailed algorithm.

- [ ] **Step 1: Write the new BrukerTimsFile.cpp**

Replace the entire file. The structure:

1. **Includes**: Replace `timsrust_cpp_bridge.h` with opentims headers + SQLiteCpp
2. **TdfMetadataReader**: Private namespace with SQL helper functions:
   - `readGlobalMetadata()` — calibration params
   - `readFrameCounts()` — MS1/MS2 frame counts
   - `readDDAStruct()` — precursor-to-frame mappings from PasefFrameMsMsInfo + Precursors
   - `readDIAWindows()` — SWATH windows from DiaFrameMsMsWindows (fallback: DiaFrameMsMsInfo)
   - `isDIA()` — detect DIA by table presence
3. **Helper structs**: Keep `ImsPeak`, `FrameCentroider`, `isCentroidingEnabled`, `expandScanOffsets` (unchanged — they operate on arrays, not timsrust types)
4. **Converter registration**: Register `OpenSourceTof2MzConverterFactory` and `OpenSourceScan2ImConverterFactory` before opening `TimsDataHandle`
5. **load()**: Open `TimsDataHandle`, detect DDA/DIA, dispatch
6. **loadDDA_()**:
   - MS1: iterate frames where `msms_type == 0`, use `save_to_buffs()` for m/z+IM+intensity
   - MS2: use `readDDAStruct()` results, extract peaks per precursor from frame scan offsets
7. **loadDIA_()**:
   - MS1: same as DDA
   - MS2: iterate frames where `msms_type != 0`, split by SWATH window IM bounds
8. **loadFrames_()**: iterate all frames, build CONCATENATED spectra
9. **transform()**: same dispatch pattern, feed to consumer

Key implementation detail for DDA MS2 reconstruction:
```cpp
// For each precursor from SQL:
//   For each frame the precursor appears in:
//     frame.save_to_buffs(frame_ids, scan_ids, tofs, intensities, mzs, ims, rts)
//     Use scan_offsets to filter peaks in [ScanNumBegin, ScanNumEnd)
//     Merge across frames, compute intensity-weighted mean IM
```

The implementation must use opentims's `TimsFrame::save_to_buffs()` which writes into caller-provided buffers. The buffers include: `frame_id`, `scan_id`, `tof`, `intensity`, `mz`, `inv_ion_mobility`, `retention_time` — all as pre-allocated arrays. The `mz` and `inv_ion_mobility` arrays are filled via the registered converters automatically.

**OLS Recalibration (when Config::calibrate = true):**
After opening the `TimsDataHandle` and loading DDA precursors from SQL, but before building MS2 spectra:
1. For each DDA precursor, get `MonoisotopicMz` from the `Precursors` table
2. In the precursor's frame, locate the highest-intensity peak within `[ScanNumBegin, ScanNumEnd)` using `scan_offsets`
3. Get that peak's raw TOF index (from the frame's `tof` buffer before conversion)
4. Convert TOF → m/z using current calibration
5. If `|converted_mz - MonoisotopicMz| < calibration_tolerance` (default 0.1 Da), add `(tof_index, sqrt(MonoisotopicMz))` to regression set
6. With >= 2 pairs, compute OLS: `slope = Cov(tof, sqrt_mz) / Var(tof)`, `intercept = mean(sqrt_mz) - slope * mean(tof)`
7. Call `converter->updateCalibration(intercept, slope)` on the `OpenSourceTof2MzConverter`
8. If < 2 pairs, log a warning and keep original calibration

This must happen before the main DDA spectrum assembly loop so all spectra use the recalibrated values.

- [ ] **Step 2: Verify compilation**

```bash
cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -20
```

Expected: compiles without errors (requires opentims to be fetchable — set the commit hash first).

- [ ] **Step 3: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat: rewrite BrukerTimsFile against opentims API with SQL-based metadata reading"
```

---

## Task 6: Update FileHandler.cpp and FileConverter.cpp

**Files:**
- Modify: `src/openms/source/FORMAT/FileHandler.cpp:898-914`
- Modify: `src/topp/FileConverter.cpp:180-241`

- [ ] **Step 1: Update FileHandler.cpp**

Replace `#ifdef WITH_TIMSRUST` with `#ifdef WITH_OPENTIMS` at three locations:
- Lines 54-56: `#include <OpenMS/FORMAT/BrukerTimsFile.h>` guard
- Lines 898 and 914: `case FileTypes::BRUKER_TDF:` block

No other changes — the `BrukerTimsFile` class name and API are unchanged.

- [ ] **Step 2: Update FileConverter.cpp TOPP parameters**

In `src/topp/FileConverter.cpp`, replace lines 180-204:

```cpp
#ifdef WITH_OPENTIMS
    registerTOPPSubsection_("bruker", "Options for reading Bruker TimsTOF .d files (requires WITH_OPENTIMS)");
    registerDoubleOption_("bruker:calibration_tolerance", "<float>", 0.0, "m/z recalibration tolerance (0 = default 0.1 Da)", false, true);
    setMinFloat_("bruker:calibration_tolerance", 0.0);
    registerStringOption_("bruker:calibrate", "<toggle>", "false", "Enable OLS m/z recalibration (may fail on some datasets)", false, true);
    setValidStrings_("bruker:calibrate", {"true", "false"});
    registerStringOption_("bruker:export_mode", "<mode>", "auto", "Export mode: 'auto' detects DDA/DIA acquisition type, "
      "'spectrum' forces per-precursor MS2 spectra (DDA-style), 'frame' returns raw 4D frames.", false, true);
    setValidStrings_("bruker:export_mode", {"auto", "spectrum", "frame"});
    registerDoubleOption_("bruker:ms1_centroid_mz_ppm", "<float>", 0.0,
      "MS1 frame IM-centroiding m/z tolerance in ppm. Both this and ms1_centroid_im_pct must be > 0 to enable. "
      "Suggested value: 5.0.", false, true);
    setMinFloat_("bruker:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:ms1_centroid_im_pct", "<float>", 0.0,
      "MS1 frame IM-centroiding ion mobility tolerance in percent. Both this and ms1_centroid_mz_ppm "
      "must be > 0 to enable. Suggested value: 3.0.", false, true);
    setMinFloat_("bruker:ms1_centroid_im_pct", 0.0);
#endif
```

Replace lines 225-241 (the `getTimsConfig_()` helper):

```cpp
#ifdef WITH_OPENTIMS
  BrukerTimsFile::Config getBrukerConfig_()
  {
    BrukerTimsFile::Config c;
    c.calibration_tolerance = getDoubleOption_("bruker:calibration_tolerance");
    c.calibrate = (getStringOption_("bruker:calibrate") == "true");
    String mode = getStringOption_("bruker:export_mode");
    if (mode == "spectrum") c.export_mode = BrukerTimsFile::Config::SPECTRUM;
    else if (mode == "frame") c.export_mode = BrukerTimsFile::Config::FRAME;
    else c.export_mode = BrukerTimsFile::Config::AUTO;
    c.ms1_centroid_mz_ppm = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_mz_ppm"));
    c.ms1_centroid_im_pct = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_im_pct"));
    return c;
  }
#endif
```

Also update all four `#ifdef WITH_TIMSRUST` locations in FileConverter.cpp:
- Line 36: `#include` guard
- Line 180: parameter registration block
- Line 225: `getTimsConfig_()` definition (rename to `getBrukerConfig_()`)
- Line 503: call site in `main_()` — rename `getTimsConfig_()` to `getBrukerConfig_()`

- [ ] **Step 3: Search for any other WITH_TIMSRUST references**

```bash
grep -rn 'WITH_TIMSRUST\|getTimsConfig_' src/topp/ src/openms/
```

Fix any remaining references.

- [ ] **Step 4: Commit**

```bash
git add src/openms/source/FORMAT/FileHandler.cpp src/topp/FileConverter.cpp
git commit -m "refactor: rename timsrust references to opentims/bruker in FileHandler and FileConverter"
```

---

## Task 7: Update Tests

**Files:**
- Modify: `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`
- Modify: `src/tests/class_tests/openms/CMakeLists.txt:154-159`
- Modify: `src/tests/topp/CMakeLists.txt`

- [ ] **Step 1: Update BrukerTimsFile_test.cpp — rename symbols**

Replace all `WITH_TIMSRUST` with `WITH_OPENTIMS`. Replace all `TIMSRUST_DDA_TEST_DATA` with `OPENTIMS_DDA_TEST_DATA`, `TIMSRUST_DIA_TEST_DATA` with `OPENTIMS_DIA_TEST_DATA`, `TIMSRUST_TEST_FASTA` with `OPENTIMS_TEST_FASTA`.

- [ ] **Step 2: Update src/tests/class_tests/openms/CMakeLists.txt**

Replace lines 154-159 — this is where the `TIMSRUST_*` compile definitions are passed to the test binary:

```cmake
# OLD:
if (TARGET BrukerTimsFile_test AND TIMSRUST_DDA_TEST_DATA)
  target_compile_definitions(BrukerTimsFile_test PRIVATE
    TIMSRUST_DDA_TEST_DATA="${TIMSRUST_DDA_TEST_DATA}"
    TIMSRUST_DIA_TEST_DATA="${TIMSRUST_DIA_TEST_DATA}"
    TIMSRUST_TEST_FASTA="${TIMSRUST_TEST_FASTA}")
endif()

# NEW:
if (TARGET BrukerTimsFile_test AND OPENTIMS_DDA_TEST_DATA)
  target_compile_definitions(BrukerTimsFile_test PRIVATE
    OPENTIMS_DDA_TEST_DATA="${OPENTIMS_DDA_TEST_DATA}"
    OPENTIMS_DIA_TEST_DATA="${OPENTIMS_DIA_TEST_DATA}"
    OPENTIMS_TEST_FASTA="${OPENTIMS_TEST_FASTA}")
endif()
```

Without this change, integration tests silently skip because the `#ifdef OPENTIMS_DDA_TEST_DATA` guards in the test file will never be defined.

- [ ] **Step 3: Add new converter and factory unit tests to BrukerTimsFile_test.cpp**

Add the following test sections (inside `#ifdef WITH_OPENTIMS` but outside the `#ifdef OPENTIMS_DDA_TEST_DATA` block, so they run even without test data):

```cpp
START_SECTION(OpenSourceTof2MzConverter unit test)
{
  // Known parameters: mz_min=100, mz_max=1700, tof_max=400000, not otofControl
  // intercept = sqrt(100) = 10.0
  // slope = (sqrt(1700) - sqrt(100)) / 400000 = (41.2310... - 10.0) / 400000
  OpenSourceTof2MzConverter conv(100.0, 1700.0, 400000, false);

  double mz_out = 0;
  uint32_t tof_in = 0;  // TOF=0 should give mz_min=100
  conv.convert(0, &mz_out, &tof_in, 1);
  TEST_REAL_SIMILAR(mz_out, 100.0);

  tof_in = 400000;  // TOF=max should give mz_max=1700
  conv.convert(0, &mz_out, &tof_in, 1);
  TEST_REAL_SIMILAR(mz_out, 1700.0);

  // Inverse: mz=100 -> tof=0
  double tof_out = 0;
  double mz_in = 100.0;
  conv.inverse_convert(0, &tof_out, &mz_in, 1);
  TEST_REAL_SIMILAR(tof_out, 0.0);
}
END_SECTION

START_SECTION(OpenSourceScan2ImConverter unit test)
{
  // Known parameters: im_min=0.6, im_max=1.6, scan_max=1000
  // intercept = 1.6 (im_max), slope = (0.6-1.6)/1000 = -0.001
  OpenSourceScan2ImConverter conv(0.6, 1.6, 1000);

  double im_out = 0;
  uint32_t scan_in = 0;  // scan=0 -> im_max=1.6
  conv.convert(0, &im_out, &scan_in, 1);
  TEST_REAL_SIMILAR(im_out, 1.6);

  scan_in = 1000;  // scan=max -> im_min=0.6
  conv.convert(0, &im_out, &scan_in, 1);
  TEST_REAL_SIMILAR(im_out, 0.6);
}
END_SECTION
```

- [ ] **Step 4: Update src/tests/topp/CMakeLists.txt**

Replace all `WITH_TIMSRUST` with `WITH_OPENTIMS`, `TIMSRUST_DDA_TEST_DATA` with `OPENTIMS_DDA_TEST_DATA`, `TIMSRUST_DIA_TEST_DATA` with `OPENTIMS_DIA_TEST_DATA`, `TIMSRUST_TEST_FASTA` with `OPENTIMS_TEST_FASTA`, `TIMSRUST_DDA_SYMLINK` with `OPENTIMS_DDA_SYMLINK`.

If any test command lines include `-timsrust:` parameters, rename them to `-bruker:`.

- [ ] **Step 5: Commit**

```bash
git add src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp \
        src/tests/class_tests/openms/CMakeLists.txt \
        src/tests/topp/CMakeLists.txt
git commit -m "test: rename timsrust to opentims in test infrastructure, add converter unit tests"
```

---

## Task 8: Build, Test, and Validate

- [ ] **Step 1: Set opentims commit hash**

After license clarification, set the actual commit hash in `cmake/cmake_findExternalLibs.cmake` (replace `PLACEHOLDER_COMMIT_HASH`).

- [ ] **Step 2: Full build**

```bash
cmake -S . -B OpenMS-build -DWITH_OPENTIMS=ON -DENABLE_OPENTIMS_TESTS=ON
cmake --build OpenMS-build -j$(nproc)
```

Expected: builds cleanly with no errors.

- [ ] **Step 3: Run unit tests**

```bash
ctest --test-dir OpenMS-build -R BrukerTimsFile -V
```

Expected: all BrukerTimsFile tests pass.

- [ ] **Step 4: Run TOPP integration tests**

```bash
ctest --test-dir OpenMS-build -R "SimpleSearchEngine_DDA|PeptideDataBaseSearchFI_DDA|SageAdapter_DDA" -V
```

Expected: all DDA TOPP tests pass.

- [ ] **Step 5: Golden-file regression**

Compare new output against golden files from Task 1. Generate new mzML from the same test data:

```bash
OpenMS-build/bin/FileConverter -in "${OPENTIMS_DDA_TEST_DATA}" \
  -out /tmp/opentims_dda.mzML
OpenMS-build/bin/FileConverter -in "${OPENTIMS_DIA_TEST_DATA}" \
  -out /tmp/opentims_dia.mzML
```

Compare spectrum counts, MS levels, precursor m/z values (allow 0.05 Th tolerance for m/z), IM values, and RT values against the golden files. Any structural differences (different spectrum counts, missing precursors) are regressions that must be fixed.

- [ ] **Step 6: Commit any fixes**

```bash
git add <specific-files-that-changed>
git commit -m "fix: address regression test findings from opentims migration"
```

---

## Task 9: Final Cleanup

- [ ] **Step 1: Verify no remaining timsrust references**

```bash
grep -rn 'timsrust\|TIMSRUST\|WITH_TIMSRUST\|tims_open\|tims_close\|tims_get_' \
  --include='*.cpp' --include='*.h' --include='*.cmake' --include='CMakeLists.txt' \
  src/ cmake/ CMakeLists.txt | grep -v 'docs/' | grep -v 'golden'
```

Expected: no matches outside of docs/specs.

- [ ] **Step 2: Verify WITH_OPENTIMS=OFF still compiles**

```bash
cmake -S . -B OpenMS-build-noopentims -DWITH_OPENTIMS=OFF
cmake --build OpenMS-build-noopentims -j$(nproc)
```

Expected: builds cleanly without opentims.

- [ ] **Step 3: Final commit**

```bash
git add <specific-files-that-changed>
git commit -m "chore: final cleanup of opentims migration"
```
