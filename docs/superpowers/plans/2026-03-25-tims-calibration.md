# Tiered TIMS Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a tiered scan-to-1/K0 calibration system (Bruker SDK → rational function → linear) to `BrukerTimsFile`, with per-frame calibration from the `TimsCalibration` table.

**Architecture:** New `RationalScan2ImConverter` class inherits from opentims's `Scan2InvIonMobilityConverter` base class. `BrukerTimsFile::Config` gains a `TimsCalibrationStrategy` enum and pressure compensation settings. The existing `openTimsDataHandle()` helper implements a tiered fallback that creates the handle with linear converters then upgrades by replacing the public `scan2inv_ion_mobility_converter` unique_ptr on `TimsDataHandle`.

**Tech Stack:** C++17, opentims (FetchContent dependency), SQLiteCpp, OpenMS test framework (ClassTest.h)

**Spec:** `docs/superpowers/specs/2026-03-25-tims-calibration-design.md`

---

## File Structure

| Action | File | Responsibility |
|--------|------|----------------|
| Create | `src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h` | Class declaration: `Coefficients` struct, converter class |
| Create | `src/openms/source/FORMAT/RationalScan2ImConverter.cpp` | Implementation: formula, inverse, `getCalibration()`, `tryCreateRationalConverter()` |
| Modify | `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h:47-57` | Add `TimsCalibrationStrategy`, `PressureCompensation` enums, `bruker_sdk_path` to Config |
| Modify | `src/openms/source/FORMAT/BrukerTimsFile.cpp:503-521` | Tiered `openTimsDataHandle()` taking Config |
| Modify | `src/openms/source/FORMAT/BrukerTimsFile.cpp:636-690` | Pass Config to `openTimsDataHandle()` at call sites |
| Modify | `src/openms/source/FORMAT/sources.cmake:126-128` | Add `RationalScan2ImConverter.cpp` |
| Modify | `src/openms/include/OpenMS/FORMAT/sources.cmake:137-139` | Add `RationalScan2ImConverter.h` |
| Modify | `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp` | Add unit tests for rational converter |

---

### Task 1: Add Config enums to BrukerTimsFile.h

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h:47-57`

- [ ] **Step 1: Add the new enums and fields to Config**

In `BrukerTimsFile.h`, replace the `Config` struct (lines 47-57) with:

```cpp
    /// Processing and export configuration
    struct Config
    {
      double calibration_tolerance = 0.0;  ///< m/z recalibration tolerance in Da (0 = default 0.1 Da)
      bool calibrate = false;              ///< Enable m/z recalibration (off by default; may fail on some datasets)

      float ms1_centroid_mz_ppm = 0.0f; ///< MS1 IM-centroiding m/z tolerance in ppm (0 = disabled, suggested: 5.0). Adapted from Sage (Lazear 2023, doi:10.1021/acs.jproteome.3c00486).
      float ms1_centroid_im_pct = 0.0f;  ///< MS1 IM-centroiding ion mobility tolerance in percent (0 = disabled, suggested: 3.0)

      enum ExportMode { AUTO, SPECTRUM, FRAME };
      ExportMode export_mode = AUTO;       ///< AUTO detects DDA vs DIA; SPECTRUM forces per-precursor; FRAME returns raw 4D frames

      /// Strategy for converting TIMS scan indices to 1/K0 values.
      /// AUTO (default): tries Bruker SDK → rational (TimsCalibration table) → linear.
      enum class TimsCalibrationStrategy { AUTO, BRUKER_SDK, RATIONAL, LINEAR };
      TimsCalibrationStrategy tims_calibration_strategy = TimsCalibrationStrategy::AUTO;

      /// Pressure compensation strategy (only effective with BRUKER_SDK calibration).
      /// Corrects for ambient gas pressure drift during acquisition.
      /// Ignored when using RATIONAL or LINEAR strategies; a warning is logged.
      enum class PressureCompensation { NONE, GLOBAL, PER_FRAME };
      PressureCompensation pressure_compensation = PressureCompensation::NONE;

      /// Path to Bruker SDK library (timsdata.dll / libtimsdata.so).
      /// Empty string (default): discover from OPENMS_BRUKER_SDK_PATH env var.
      std::string bruker_sdk_path;
    };
```

Add `#include <string>` to the includes if not already present.

- [ ] **Step 2: Verify build**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -20`
Expected: builds cleanly (Config is header-only, no ABI change for default-constructed configs)

- [ ] **Step 3: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/BrukerTimsFile.h
git commit -m "feat(BrukerTimsFile): add TimsCalibrationStrategy and PressureCompensation to Config"
```

---

### Task 2: Create RationalScan2ImConverter header

**Files:**
- Create: `src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h`
- Modify: `src/openms/include/OpenMS/FORMAT/sources.cmake:137-139`

- [ ] **Step 1: Write the header**

Create `src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/OpenMSConfig.h>
#include <opentims++/scan2inv_ion_mobility_converter.h>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Per-frame rational function converter using Bruker TimsCalibration table.
   *
   * This is the first open-source implementation to read per-frame calibration
   * coefficients from the TimsCalibration table in analysis.tdf. All other
   * open-source tools (OpenTIMS, AlphaTIMS, timsrust) use a simpler linear
   * approximation from GlobalMetadata that is frame-independent.
   *
   * Physical model (ModelType=2):
   *   Inner term maps scan index to TIMS elution voltage:
   *     V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)
   *   Outer term maps voltage to inverse reduced ion mobility:
   *     1/K0 = 1 / (c6 + c7 / V)
   *
   * Combined:
   *   1/K0 = 1 / (c6 + c7 / (c2 + ((c3 - c2) / c1) * (scan - c4 - c0)))
   *
   * Limitations vs Bruker SDK:
   *   - No pressure compensation (ambient gas pressure drift correction)
   *   - No recalibrated-state awareness (uses raw calibration from acquisition)
   *   - Only ModelType=2 supported; unknown types should fall back to linear
   *
   * Thread safety: immutable after construction. convert() and inverse_convert()
   * are safe for concurrent calls.
   */
  class OPENMS_DLLAPI RationalScan2ImConverter : public Scan2InvIonMobilityConverter
  {
  public:
    /// Calibration coefficients from one row of the TimsCalibration table.
    struct Coefficients
    {
      double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9;
      // c5, c8, c9: stored but unused in ModelType=2. No public documentation
      // assigns them meaning. Retained for forward compatibility.
    };

    /// @param calibrations  Map of calibration ID -> coefficients (from TimsCalibration table)
    /// @param frame_to_cal  Calibration ID for each frame, indexed by frame_id (1-based; index 0 unused)
    RationalScan2ImConverter(
      std::unordered_map<uint32_t, Coefficients> calibrations,
      std::vector<uint32_t> frame_to_cal);

    void convert(uint32_t frame_id, double* inv_ion_mobilities,
                 const double* scans, uint32_t size) override;
    void convert(uint32_t frame_id, double* inv_ion_mobilities,
                 const uint32_t* scans, uint32_t size) override;
    void inverse_convert(uint32_t frame_id, uint32_t* scans,
                         const double* inv_ion_mobilities, uint32_t size) override;

    /// Returns e.g. "RationalScan2ImConverter (2 calibration segments)"
    std::string description() const override;

  private:
    std::unordered_map<uint32_t, Coefficients> calibrations_;
    std::vector<uint32_t> frame_to_cal_;  ///< indexed by frame_id (1-based)

    /// Look up calibration for a frame. Falls back to first calibration for
    /// out-of-range frame_id with a warning.
    const Coefficients& getCalibration(uint32_t frame_id) const;

    /// V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)
    /// 1/K0 = 1.0 / (c6 + c7 / V)
    static double applyFormula(const Coefficients& c, double scan);

    /// scan = c0 + c4 + (c1 / (c3 - c2)) * (c7 / (1.0/inv_k0 - c6) - c2)
    static double invertFormula(const Coefficients& c, double inv_k0);
  };

  /// Try to create a RationalScan2ImConverter by reading the TimsCalibration
  /// table from the given .d directory. Returns nullptr if the table is missing,
  /// has unknown ModelType, or any frame has a NULL/invalid TimsCalibration FK.
  OPENMS_DLLAPI std::unique_ptr<Scan2InvIonMobilityConverter> tryCreateRationalConverter(
    const std::string& tims_dir_path);

} // namespace OpenMS

#endif // WITH_OPENTIMS
```

- [ ] **Step 2: Register header in sources.cmake**

In `src/openms/include/OpenMS/FORMAT/sources.cmake`, after line 138 (`list(APPEND sources_list_h BrukerTimsFile.h)`), add:

```cmake
  list(APPEND sources_list_h RationalScan2ImConverter.h)
```

- [ ] **Step 3: Verify build**

Run: `cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -20`
Expected: builds (header-only so far, no .cpp yet — cmake reconfigure may be needed)

- [ ] **Step 4: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h \
        src/openms/include/OpenMS/FORMAT/sources.cmake
git commit -m "feat: add RationalScan2ImConverter header (per-frame TIMS calibration)"
```

---

### Task 3: Write unit tests for RationalScan2ImConverter

**Files:**
- Modify: `src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp`

- [ ] **Step 1: Write failing tests**

Add these tests after the existing `END_SECTION` on line 31 (after the "invalid path" test) and before the `#ifdef OPENTIMS_DDA_TEST_DATA` block on line 46. These tests do not require test data:

```cpp
START_SECTION(RationalScan2ImConverter forward conversion)
{
  // Coefficients from opentims test data (ModelType=2)
  using Coeff = RationalScan2ImConverter::Coefficients;
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  std::unordered_map<uint32_t, Coeff> calibrations;
  calibrations[1] = c;
  std::vector<uint32_t> frame_to_cal = {0, 1}; // index 0 unused, frame 1 -> cal 1

  RationalScan2ImConverter converter(std::move(calibrations), std::move(frame_to_cal));

  // Worked example from spec (precise values):
  // V = 213.5998 + ((75.81729 - 213.5998) / 917.0) * (500 - 33 - 1) = 143.5816
  // 1/K0 = 1 / (-0.009065829 + 135.4364 / 143.5816) = 1.070429
  double scan_val = 500.0;
  double result = 0.0;
  converter.convert(1, &result, &scan_val, 1);
  TEST_REAL_SIMILAR(result, 1.070429);

  // Test with uint32_t scan input
  uint32_t scan_int = 500;
  double result2 = 0.0;
  converter.convert(1, &result2, &scan_int, 1);
  TEST_REAL_SIMILAR(result2, 1.070429);

  // Test multiple scans at once
  double scans[] = {0.0, 250.0, 500.0, 916.0};
  double results[4] = {};
  converter.convert(1, results, scans, 4);
  // scan=0: V = 213.5998 + slope*(0-34) = 218.7084
  //         1/K0 = 1/(-0.009066 + 135.4364/218.7084) = 1.63883
  TEST_REAL_SIMILAR(results[0], 1.63883);
  // All should be positive
  for (int i = 0; i < 4; ++i)
  {
    TEST_EQUAL(results[i] > 0.0, true);
  }
  // Increasing scan -> decreasing 1/K0 (higher scan = lower mobility)
  for (int i = 1; i < 4; ++i)
  {
    TEST_EQUAL(results[i] < results[i-1], true);
  }
}
END_SECTION

START_SECTION(RationalScan2ImConverter round-trip via inverse_convert)
{
  using Coeff = RationalScan2ImConverter::Coefficients;
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  std::unordered_map<uint32_t, Coeff> calibrations;
  calibrations[1] = c;
  std::vector<uint32_t> frame_to_cal = {0, 1};

  RationalScan2ImConverter converter(std::move(calibrations), std::move(frame_to_cal));

  // Forward: scan 500 -> 1/K0
  double scan_val = 500.0;
  double inv_k0 = 0.0;
  converter.convert(1, &inv_k0, &scan_val, 1);

  // Inverse: 1/K0 -> scan (should round-trip to 500)
  uint32_t scan_back = 0;
  converter.inverse_convert(1, &scan_back, &inv_k0, 1);
  TEST_EQUAL(scan_back, 500);

  // Test a range of scans for round-trip
  for (uint32_t s = 10; s < 900; s += 50)
  {
    double sv = static_cast<double>(s);
    double ik0 = 0.0;
    converter.convert(1, &ik0, &sv, 1);
    uint32_t back = 0;
    converter.inverse_convert(1, &back, &ik0, 1);
    // Allow +/- 1 for rounding
    TEST_EQUAL(back >= s - 1 && back <= s + 1, true);
  }
}
END_SECTION

START_SECTION(RationalScan2ImConverter per-frame calibration)
{
  using Coeff = RationalScan2ImConverter::Coefficients;
  // Two different calibration segments
  Coeff c1{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};
  Coeff c2{1.0, 917.0, 220.0, 80.0, 33.0, 1.0, -0.01, 140.0, 13.0, 1660.0};

  std::unordered_map<uint32_t, Coeff> calibrations;
  calibrations[1] = c1;
  calibrations[2] = c2;
  // frame 1 uses cal 1, frame 2 uses cal 2
  std::vector<uint32_t> frame_to_cal = {0, 1, 2};

  RationalScan2ImConverter converter(std::move(calibrations), std::move(frame_to_cal));

  double scan_val = 500.0;
  double result_f1 = 0.0, result_f2 = 0.0;
  converter.convert(1, &result_f1, &scan_val, 1);
  converter.convert(2, &result_f2, &scan_val, 1);

  // Different calibrations should produce different results
  TEST_NOT_EQUAL(result_f1, result_f2);
  // Both should be positive
  TEST_EQUAL(result_f1 > 0.0, true);
  TEST_EQUAL(result_f2 > 0.0, true);
}
END_SECTION

START_SECTION(RationalScan2ImConverter description)
{
  using Coeff = RationalScan2ImConverter::Coefficients;
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  std::unordered_map<uint32_t, Coeff> cals;
  cals[1] = c;
  std::vector<uint32_t> ftc = {0, 1};

  RationalScan2ImConverter converter(std::move(cals), std::move(ftc));
  std::string desc = converter.description();
  TEST_EQUAL(desc.find("RationalScan2ImConverter") != std::string::npos, true);
  TEST_EQUAL(desc.find("1") != std::string::npos, true); // 1 calibration segment
}
END_SECTION

START_SECTION(RationalScan2ImConverter singularity edge cases)
{
  using Coeff = RationalScan2ImConverter::Coefficients;

  // Normal coefficients as baseline
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};

  // Edge case: c3 == c2 (zero slope in scan-to-voltage mapping)
  // Should not crash — singularity guard produces a finite value
  Coeff c_zero_slope = c;
  c_zero_slope.c3 = c_zero_slope.c2; // c3 == c2

  std::unordered_map<uint32_t, Coeff> cals;
  cals[1] = c_zero_slope;
  std::vector<uint32_t> ftc = {0, 1};
  RationalScan2ImConverter conv(std::move(cals), std::move(ftc));

  double scan_val = 500.0;
  double result = 0.0;
  conv.convert(1, &result, &scan_val, 1);
  TEST_EQUAL(std::isfinite(result), true);

  // Edge case: scan that makes V approach 0
  // With normal coefficients, V=0 is physically unreachable, but we test the guard
  // by using coefficients where V can be 0 for some scan value
  Coeff c_v_zero = c;
  c_v_zero.c2 = 0.0; // offset = 0
  c_v_zero.c3 = 0.0; // slope endpoint = 0 (V always 0 regardless of scan)

  std::unordered_map<uint32_t, Coeff> cals2;
  cals2[1] = c_v_zero;
  std::vector<uint32_t> ftc2 = {0, 1};
  RationalScan2ImConverter conv2(std::move(cals2), std::move(ftc2));

  double result2 = 0.0;
  conv2.convert(1, &result2, &scan_val, 1);
  TEST_EQUAL(std::isfinite(result2), true);
}
END_SECTION
```

Also add the include at the top of the file (after line 12):

```cpp
#include <OpenMS/FORMAT/RationalScan2ImConverter.h>
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc) 2>&1 | tail -20`
Expected: FAIL — `RationalScan2ImConverter` methods not defined yet (linker errors)

- [ ] **Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp
git commit -m "test: add unit tests for RationalScan2ImConverter (red phase)"
```

---

### Task 4: Implement RationalScan2ImConverter

**Files:**
- Create: `src/openms/source/FORMAT/RationalScan2ImConverter.cpp`
- Modify: `src/openms/source/FORMAT/sources.cmake:126-128`

- [ ] **Step 1: Write the implementation**

Create `src/openms/source/FORMAT/RationalScan2ImConverter.cpp`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/RationalScan2ImConverter.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <SQLiteCpp/SQLiteCpp.h>
#include <cmath>
#include <sstream>

namespace OpenMS
{

  RationalScan2ImConverter::RationalScan2ImConverter(
    std::unordered_map<uint32_t, Coefficients> calibrations,
    std::vector<uint32_t> frame_to_cal)
    : calibrations_(std::move(calibrations)),
      frame_to_cal_(std::move(frame_to_cal))
  {
  }

  const RationalScan2ImConverter::Coefficients&
  RationalScan2ImConverter::getCalibration(uint32_t frame_id) const
  {
    // frame_to_cal_ is 1-based (index 0 unused)
    if (frame_id > 0 && frame_id < frame_to_cal_.size())
    {
      auto it = calibrations_.find(frame_to_cal_[frame_id]);
      if (it != calibrations_.end())
      {
        return it->second;
      }
    }
    // Fallback: use first calibration entry (should not happen in valid data)
    if (frame_id != 0)
    {
      OPENMS_LOG_WARN << "RationalScan2ImConverter: frame_id " << frame_id
                      << " out of range, using first calibration" << std::endl;
    }
    return calibrations_.begin()->second;
  }

  double RationalScan2ImConverter::applyFormula(const Coefficients& c, double scan)
  {
    // V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)
    double V = c.c2 + ((c.c3 - c.c2) / c.c1) * (scan - c.c4 - c.c0);
    if (V == 0.0) V = 1e-10; // guard against division by zero (physically unreachable)
    // 1/K0 = 1 / (c6 + c7 / V)
    double denom = c.c6 + c.c7 / V;
    if (denom == 0.0) denom = 1e-10; // guard
    return 1.0 / denom;
  }

  double RationalScan2ImConverter::invertFormula(const Coefficients& c, double inv_k0)
  {
    // scan = c0 + c4 + (c1 / (c3 - c2)) * (c7 / (1/inv_k0 - c6) - c2)
    double recip = 1.0 / inv_k0;
    double inner_denom = recip - c.c6;
    if (inner_denom == 0.0) inner_denom = 1e-10; // guard
    double slope_denom = c.c3 - c.c2;
    if (slope_denom == 0.0) slope_denom = 1e-10; // guard
    return c.c0 + c.c4 + (c.c1 / slope_denom) * (c.c7 / inner_denom - c.c2);
  }

  void RationalScan2ImConverter::convert(uint32_t frame_id,
    double* inv_ion_mobilities, const double* scans, uint32_t size)
  {
    const auto& c = getCalibration(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = applyFormula(c, scans[i]);
    }
  }

  void RationalScan2ImConverter::convert(uint32_t frame_id,
    double* inv_ion_mobilities, const uint32_t* scans, uint32_t size)
  {
    const auto& c = getCalibration(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = applyFormula(c, static_cast<double>(scans[i]));
    }
  }

  void RationalScan2ImConverter::inverse_convert(uint32_t frame_id,
    uint32_t* scans, const double* inv_ion_mobilities, uint32_t size)
  {
    const auto& c = getCalibration(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      double val = invertFormula(c, inv_ion_mobilities[i]);
      scans[i] = val > 0.0 ? static_cast<uint32_t>(val + 0.5) : 0;
    }
  }

  std::string RationalScan2ImConverter::description() const
  {
    return "RationalScan2ImConverter (" + std::to_string(calibrations_.size())
           + " calibration segment" + (calibrations_.size() != 1 ? "s" : "") + ")";
  }

  // =========================================================================
  // Factory function: try to build a RationalScan2ImConverter from SQLite
  // =========================================================================

  std::unique_ptr<Scan2InvIonMobilityConverter> tryCreateRationalConverter(
    const std::string& tims_dir_path)
  {
    std::string tdf_path = tims_dir_path + "/analysis.tdf";

    try
    {
      SQLite::Database db(tdf_path, SQLite::OPEN_READONLY);

      // 1. Read TimsCalibration table
      std::unordered_map<uint32_t, RationalScan2ImConverter::Coefficients> calibrations;
      {
        SQLite::Statement q(db,
          "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 "
          "FROM TimsCalibration");

        while (q.executeStep())
        {
          uint32_t id = static_cast<uint32_t>(q.getColumn(0).getInt());
          int model_type = q.getColumn(1).getInt();

          if (model_type != 2)
          {
            OPENMS_LOG_WARN << "TimsCalibration ModelType=" << model_type
                            << " unsupported (only ModelType=2 implemented), "
                            << "falling back to linear converter" << std::endl;
            return nullptr;
          }

          RationalScan2ImConverter::Coefficients c;
          c.c0 = q.getColumn(2).getDouble();
          c.c1 = q.getColumn(3).getDouble();
          c.c2 = q.getColumn(4).getDouble();
          c.c3 = q.getColumn(5).getDouble();
          c.c4 = q.getColumn(6).getDouble();
          c.c5 = q.getColumn(7).getDouble();
          c.c6 = q.getColumn(8).getDouble();
          c.c7 = q.getColumn(9).getDouble();
          c.c8 = q.getColumn(10).getDouble();
          c.c9 = q.getColumn(11).getDouble();

          calibrations[id] = c;
        }
      }

      if (calibrations.empty())
      {
        OPENMS_LOG_DEBUG << "TimsCalibration table empty, "
                         << "falling back to linear converter" << std::endl;
        return nullptr;
      }

      // 2. Read frame-to-calibration mapping from Frames table
      //    (separate try/catch: column may not exist in old TDF versions)
      std::vector<uint32_t> frame_to_cal;
      try
      {
        SQLite::Statement q(db, "SELECT Id, TimsCalibration FROM Frames ORDER BY Id");
        while (q.executeStep())
        {
          uint32_t frame_id = static_cast<uint32_t>(q.getColumn(0).getInt());
          if (q.getColumn(1).isNull())
          {
            OPENMS_LOG_WARN << "Frame " << frame_id
                            << " has NULL TimsCalibration FK, "
                            << "falling back to linear converter" << std::endl;
            return nullptr;
          }
          uint32_t cal_id = static_cast<uint32_t>(q.getColumn(1).getInt());

          // Ensure vector is large enough (frame IDs are 1-based)
          if (frame_to_cal.size() <= frame_id)
          {
            frame_to_cal.resize(frame_id + 1, 0);
          }
          frame_to_cal[frame_id] = cal_id;

          // Verify calibration ID exists
          if (calibrations.find(cal_id) == calibrations.end())
          {
            OPENMS_LOG_WARN << "Frame " << frame_id
                            << " references unknown TimsCalibration Id=" << cal_id
                            << ", falling back to linear converter" << std::endl;
            return nullptr;
          }
        }
      }
      catch (const SQLite::Exception&)
      {
        // TimsCalibration column does not exist in Frames table (old TDF format)
        OPENMS_LOG_DEBUG << "Frames table has no TimsCalibration column, "
                         << "falling back to linear converter" << std::endl;
        return nullptr;
      }

      if (frame_to_cal.size() <= 1) // only index 0 (unused)
      {
        OPENMS_LOG_DEBUG << "No frames found, falling back to linear converter"
                         << std::endl;
        return nullptr;
      }

      OPENMS_LOG_INFO << "TIMS calibration: rational (TimsCalibration table, "
                      << calibrations.size() << " segment"
                      << (calibrations.size() != 1 ? "s" : "")
                      << " for " << (frame_to_cal.size() - 1) << " frames)"
                      << std::endl;

      return std::make_unique<RationalScan2ImConverter>(
        std::move(calibrations), std::move(frame_to_cal));
    }
    catch (const SQLite::Exception&)
    {
      // TimsCalibration table does not exist
      OPENMS_LOG_DEBUG << "TimsCalibration table not found in " << tdf_path
                       << ", falling back to linear converter" << std::endl;
      return nullptr;
    }
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
```

- [ ] **Step 2: Register source in sources.cmake**

In `src/openms/source/FORMAT/sources.cmake`, after line 127 (`list(APPEND sources_list BrukerTimsFile.cpp)`), add:

```cmake
  list(APPEND sources_list RationalScan2ImConverter.cpp)
```

- [ ] **Step 3: Build and run tests**

Run cmake reconfigure then build:
```bash
cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc) 2>&1 | tail -20
```
Expected: builds successfully

Run tests:
```bash
ctest --test-dir OpenMS-build -R BrukerTimsFile_test -V 2>&1 | tail -40
```
Expected: all new RationalScan2ImConverter tests PASS

- [ ] **Step 4: Commit**

```bash
git add src/openms/source/FORMAT/RationalScan2ImConverter.cpp \
        src/openms/source/FORMAT/sources.cmake
git commit -m "feat: implement RationalScan2ImConverter with per-frame calibration"
```

---

### Task 5: Wire tiered fallback into BrukerTimsFile.cpp

**Files:**
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp:503-521` (openTimsDataHandle)
- Modify: `src/openms/source/FORMAT/BrukerTimsFile.cpp:636-690` (call sites)

- [ ] **Step 1: Add include for the new converter**

At the top of `BrukerTimsFile.cpp`, after line 19 (`#include <opentims++/scan2inv_ion_mobility_converter.h>`), add:

```cpp
#include <OpenMS/FORMAT/RationalScan2ImConverter.h>
```

Note: `BrukerScan2InvIonMobilityConverterFactory` is already available via the
existing `#include <opentims++/scan2inv_ion_mobility_converter.h>` on line 19.
Do NOT add `#include <opentims++/converters.h>` — it only contains
`setup_bruker()`/`setup_opensource()` which we don't use.

- [ ] **Step 2: Add pressure compensation mapping helper**

Before `openTimsDataHandle()` (before line 500), add:

```cpp
  // Map OpenMS PressureCompensation enum to opentims pressure_compensation_strategy
  static pressure_compensation_strategy mapPressureCompensation(
    BrukerTimsFile::Config::PressureCompensation pc)
  {
    switch (pc)
    {
      case BrukerTimsFile::Config::PressureCompensation::GLOBAL:
        return AnalyisGlobalPressureCompensation;
      case BrukerTimsFile::Config::PressureCompensation::PER_FRAME:
        return PerFramePressureCompensationWithMissingReference;
      default:
        return NoPressureCompensation;
    }
  }
```

- [ ] **Step 3: Rewrite openTimsDataHandle() with tiered fallback**

Replace the existing `openTimsDataHandle()` function (lines 503-521) with:

```cpp
  using Config = BrukerTimsFile::Config;

  static std::unique_ptr<TimsDataHandle> openTimsDataHandle(
    const String& path, const Config& config = Config())
  {
    std::string path_string = path;

    // 1. Always create handle with linear converters (safe baseline)
    auto& tof_factory = OpenSourceTof2MzConverterFactory::instance();
    auto& im_factory = OpenSourceScan2ImConverterFactory::instance();

    std::unique_ptr<TimsDataHandle> handle;
    try
    {
      handle = std::make_unique<TimsDataHandle>(
        path_string, NoPressureCompensation, &tof_factory, &im_factory);
    }
    catch (const std::exception& e)
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        path + " (opentims: " + String(e.what()) + ")");
    }

    using Strategy = Config::TimsCalibrationStrategy;
    auto strategy = config.tims_calibration_strategy;

    // Warn if pressure compensation requested without SDK strategy
    if (config.pressure_compensation != Config::PressureCompensation::NONE
        && strategy != Strategy::BRUKER_SDK
        && strategy != Strategy::AUTO)
    {
      OPENMS_LOG_WARN << "Pressure compensation requires BRUKER_SDK strategy, ignoring"
                      << std::endl;
    }

    // 2. Try Bruker SDK (AUTO or BRUKER_SDK)
    if (strategy == Strategy::AUTO || strategy == Strategy::BRUKER_SDK)
    {
      std::string sdk_path = config.bruker_sdk_path;
      if (sdk_path.empty())
      {
        const char* env = std::getenv("OPENMS_BRUKER_SDK_PATH");
        if (env) sdk_path = env;
      }

      if (!sdk_path.empty())
      {
        try
        {
          auto pcs = mapPressureCompensation(config.pressure_compensation);
          handle->scan2inv_ion_mobility_converter =
            BrukerScan2InvIonMobilityConverterFactory::instance(sdk_path)
              .produce(*handle, pcs);
          OPENMS_LOG_INFO << "TIMS calibration: Bruker SDK"
                          << (pcs != NoPressureCompensation ? " (pressure_comp=on)" : "")
                          << std::endl;
          return handle;
        }
        catch (const std::exception& e)
        {
          if (strategy == Strategy::BRUKER_SDK)
          {
            throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              String("Bruker SDK failed: ") + e.what());
          }
          OPENMS_LOG_DEBUG << "Bruker SDK not available (" << e.what()
                          << "), trying rational" << std::endl;
        }
      }
      else if (strategy == Strategy::BRUKER_SDK)
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Bruker SDK path not set (use OPENMS_BRUKER_SDK_PATH or Config::bruker_sdk_path)");
      }
    }

    // 3. Try rational model from TimsCalibration table (AUTO or RATIONAL)
    if (strategy == Strategy::AUTO || strategy == Strategy::RATIONAL)
    {
      auto converter = tryCreateRationalConverter(path_string);
      if (converter)
      {
        handle->scan2inv_ion_mobility_converter = std::move(converter);
        return handle;
      }
      else if (strategy == Strategy::RATIONAL)
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "TimsCalibration table not found or unsupported in " + path);
      }
    }

    // 4. Linear fallback (already set from step 1)
    OPENMS_LOG_INFO << "TIMS calibration: linear (GlobalMetadata)" << std::endl;
    return handle;
  }
```

- [ ] **Step 4: Update call sites to pass Config**

Line 644 — `load()` already has `config`, change from:
```cpp
    auto handle = openTimsDataHandle(path);
```
to:
```cpp
    auto handle = openTimsDataHandle(path, config);
```

Line 690 — `transform()` already has `config`, change from:
```cpp
    auto handle = openTimsDataHandle(path);
```
to:
```cpp
    auto handle = openTimsDataHandle(path, config);
```

The no-Config overloads (`load(path, exp)` on line 636 and `transform(path, consumer)` on line 683) already delegate to the Config overloads with `Config()`, so they are unaffected.

- [ ] **Step 5: Build and run tests**

```bash
cmake --build OpenMS-build --target BrukerTimsFile_test -j$(nproc) 2>&1 | tail -20
ctest --test-dir OpenMS-build -R BrukerTimsFile_test -V 2>&1 | tail -40
```
Expected: all tests PASS (existing + new)

- [ ] **Step 6: Commit**

```bash
git add src/openms/source/FORMAT/BrukerTimsFile.cpp
git commit -m "feat(BrukerTimsFile): wire tiered TIMS calibration fallback (SDK > rational > linear)"
```

---

### Task 6: Verify full build and existing tests

**Files:** None (verification only)

- [ ] **Step 1: Full library build**

```bash
cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -30
```
Expected: clean build, no warnings in new files

- [ ] **Step 2: Run BrukerTimsFile tests**

```bash
ctest --test-dir OpenMS-build -R BrukerTimsFile_test -V
```
Expected: all tests PASS

- [ ] **Step 3: Run broader FORMAT tests to check for regressions**

```bash
ctest --test-dir OpenMS-build -R "FORMAT" --timeout 120 2>&1 | tail -30
```
Expected: no regressions

- [ ] **Step 4: Commit (if any fixups needed)**

Only if previous steps revealed issues that required fixes.
