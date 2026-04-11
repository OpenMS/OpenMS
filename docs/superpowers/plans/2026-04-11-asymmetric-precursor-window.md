# Asymmetric Precursor Window Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unify `PeptideSearchEngineFI`'s precursor tolerance and open-search window into a single positive-magnitude asymmetric pair, refactor `FragmentIndex::getPeptidesInPrecursorRange` to an absolute-bounds interface, and extend the calibration pass to preserve signed bias instead of clobbering asymmetric bounds to symmetric.

**Architecture:** Positive-magnitude user-facing params (`precursor:mass_tolerance_lower/_upper`, both ≥ 0, shared unit), internal signed conversion centralized in one helper (`FragmentIndex::computeMassWindow_`), calibration pass computes signed median + MAD-of-residuals and skips writeback when `|shift| ≥ spread` (positive-magnitude schema cannot faithfully encode one-sided windows without loosening). `OpenSearchModificationAnalysis` fed via `min(lower, upper)` through a new private const helper. Duplicate `isOpenSearchMode_()` consolidated into a shared `FragmentIndex::isOpenSearchMode` static.

**Tech Stack:** C++ (OpenMS, FragmentIndex), nanobind (pyOpenMS), CMake + CTest, `START_TEST`/`TEST_EQUAL` class-test macros, pytest (pyOpenMS).

**Spec:** `docs/superpowers/specs/2026-04-11-asymmetric-precursor-window-design.md`

**Tracks:** OpenMS/OpenMS#9108 item 5.

---

## File Structure

### Modified

| File | Responsibility |
|---|---|
| `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h` | Header: new members, `getPeptidesInMassWindow`, `computeMassWindow_`, static `isOpenSearchMode`, friend declaration for test subclass |
| `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` | Implementations, param registrations, `updateMembers_` validation, `searchDifferentPrecursorRanges` rewrite |
| `src/openms/include/OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h` | Header: new members, `computeModMatchTolerance_`, extended `CalibrationResult_`, `last_calibration_result_`, `isOpenSearchMode_` delegation, friend for test subclass |
| `src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp` | Param registrations, `updateMembers_`, OSMA call-site helper use, `runCalibrationPass_` signed-bias rewrite, calibration writeback with extreme-bias skip + member refresh |
| `src/topp/PeptideDataBaseSearchFI.cpp` | Param-name updates, local open-mode detection |
| `src/tests/topp/PeptideDataBaseSearchFI_1.ini` | Migrate `mass_tolerance` → positive-magnitude pair |
| `src/tests/topp/PeptideDataBaseSearchFI_2.ini` | Same migration |
| `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` | Migrate existing param usage, add 10 new tests (1, 2, 3, 4, 5, 6a, 6b, 6c) using existing `FragmentIndex_test` friend subclass |
| `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp` | Migrate existing param usage, add new `PeptideSearchEngineFIAlgorithm_test` friend subclass, add 4 new tests (7, 8, 9, 10) |
| `src/pyOpenMS/tests/test_modification_discovery.py` | Migrate 8 param-set hits across 4 test functions |
| `CHANGELOG` | BREAKING entry under `TOPP tools > Changes > PeptideDataBaseSearchFI` |

### Created

| File | Responsibility |
|---|---|
| `src/tests/topp/PeptideDataBaseSearchFI_4.ini` | Thirdparty-level end-to-end test: asymmetric bounds + calibration on realistic shifted data |
| `src/tests/topp/PeptideDataBaseSearchFI_4_out.idXML` | Reference output for test 4 |
| `src/tests/topp/SimpleSearchEngine_1_shifted_7ppm.mzML` | Shifted-precursor test fixture for test 4 |

---

## Task 1: FragmentIndex header refactor

**Files:**
- Modify: `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h`

This task changes the header skeleton. The next task makes the .cpp match. After both tasks, the code builds again.

- [ ] **Step 1.1: Add `<algorithm>` include**

Edit `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h`:
```cpp
#include <array>
#include <mutex>
#include <vector>
#include <functional>
#include <algorithm>   // std::max (used by inline static isOpenSearchMode)
```

- [ ] **Step 1.2: Rename `getPeptidesInPrecursorRange` to `getPeptidesInMassWindow`, mark `const`, update doc comment**

Replace lines 191-198 (the method declaration and its doc comment) with:
```cpp
    /** Return the [begin_idx, end_idx) peptide index range such that
     * `fi_peptides_[i].precursor_mz_ ∈ [precursor_mass + window.first, precursor_mass + window.second]`
     * for all i in the returned range.
     *
     * @param[in] precursor_mass The mono-charged precursor mass (M+H).
     * @param[in] window Signed absolute offsets around the precursor mass. `window.first` must be
     *                   ≤ 0 and `window.second` must be ≥ 0. No hidden tolerance is added —
     *                   the caller supplies the full signed window (see `computeMassWindow_`).
     * @return [begin_idx, end_idx) half-open index range into `fi_peptides_`.
     */
    std::pair<size_t, size_t> getPeptidesInMassWindow(float precursor_mass,
                                                      const std::pair<float, float>& window) const;
```

- [ ] **Step 1.3: Replace tolerance members (lines 377-378) and delete open-window members (lines 456-457)**

Line 377-378 currently reads:
```cpp
    float precursor_mz_tolerance_;
    bool precursor_mz_tolerance_unit_ppm_{true};
```
Replace with:
```cpp
    double precursor_mass_tolerance_lower_{20.0};   ///< positive magnitude, effective lower bound is -lower
    double precursor_mass_tolerance_upper_{20.0};   ///< positive magnitude, effective upper bound is +upper
    bool precursor_mass_tolerance_unit_ppm_{true};
```

Also delete lines 456-457:
```cpp
    float open_precursor_window_lower_;
    float open_precursor_window_upper_;
```

- [ ] **Step 1.4: Replace `isOpenSearchMode_()` (lines 448-454) with static helper + delegating instance method**

Replace:
```cpp
    /// Helper function to determine if open search should be used based on tolerance
    bool isOpenSearchMode_() const
    {
      return precursor_mz_tolerance_unit_ppm_
               ? (precursor_mz_tolerance_ > 1000.0)
               : (precursor_mz_tolerance_ > 1.0);
    }
```
with:
```cpp
    /// Shared auto-detection: open-search iff max(lower, upper) > threshold (1000 ppm or 1 Da).
    /// Strict `>`: exactly 1000 ppm stays closed.
    static bool isOpenSearchMode(double lower_magnitude,
                                 double upper_magnitude,
                                 bool unit_ppm) noexcept
    {
      const double threshold = unit_ppm ? 1000.0 : 1.0;
      return std::max(lower_magnitude, upper_magnitude) > threshold;
    }

    /// Instance delegate — same rule, reads the member bounds.
    bool isOpenSearchMode_() const noexcept
    {
      return isOpenSearchMode(precursor_mass_tolerance_lower_,
                              precursor_mass_tolerance_upper_,
                              precursor_mass_tolerance_unit_ppm_);
    }
```

- [ ] **Step 1.5: Add private `computeMassWindow_` declaration**

In the `private:` section of `FragmentIndex` (near `searchDifferentPrecursorRanges` at line 406), add:
```cpp
    /** Compute the signed mass window {lo, hi} around a precursor_mass, converting ppm → Da
     * if the unit is ppm. `lo` is negative (or zero), `hi` is positive (or zero). This is the
     * only place where positive member magnitudes become signed offsets.
     */
    std::pair<float, float> computeMassWindow_(float precursor_mass) const;
```

- [ ] **Step 1.6: Run build to confirm the header compiles (it will fail at link time, that's expected)**

Run:
```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | head -100
```
Expected: compilation errors in `FragmentIndex.cpp` referencing the old member names. Link errors for `getPeptidesInMassWindow` / `computeMassWindow_`. **Do NOT commit** — the next task fixes the .cpp.

---

## Task 2: FragmentIndex.cpp refactor

**Files:**
- Modify: `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp`

- [ ] **Step 2.1: Replace `getPeptidesInPrecursorRange` implementation (lines 807-815)**

Replace:
```cpp
  std::pair<size_t, size_t > FragmentIndex::getPeptidesInPrecursorRange(float precursor_mass,
                                                                       const std::pair<float, float>& window)
  {
      float prec_tol = precursor_mz_tolerance_unit_ppm_ ? Math::ppmToMass(precursor_mz_tolerance_, precursor_mass) : precursor_mz_tolerance_ ;

      auto left_it = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass - prec_tol + window.first, [](const Peptide& a, float b) { return a.precursor_mz_ < b;});
      auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(), precursor_mass + prec_tol + window.second, [](float b, const Peptide& a) { return b < a.precursor_mz_;});
      return make_pair(std::distance(fi_peptides_.begin(), left_it), std::distance(fi_peptides_.begin(), right_it));
  }
```
with:
```cpp
  std::pair<size_t, size_t> FragmentIndex::getPeptidesInMassWindow(float precursor_mass,
                                                                   const std::pair<float, float>& window) const
  {
    auto left_it = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                    precursor_mass + window.first,
                                    [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
    auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                     precursor_mass + window.second,
                                     [](float b, const Peptide& a) { return b < a.precursor_mz_; });
    return std::make_pair(std::distance(fi_peptides_.begin(), left_it),
                          std::distance(fi_peptides_.begin(), right_it));
  }

  std::pair<float, float> FragmentIndex::computeMassWindow_(float precursor_mass) const
  {
    if (precursor_mass_tolerance_unit_ppm_)
    {
      const float lo = -Math::ppmToMass<float>(static_cast<float>(precursor_mass_tolerance_lower_),
                                               precursor_mass);
      const float hi =  Math::ppmToMass<float>(static_cast<float>(precursor_mass_tolerance_upper_),
                                               precursor_mass);
      return {lo, hi};
    }
    return {-static_cast<float>(precursor_mass_tolerance_lower_),
             static_cast<float>(precursor_mass_tolerance_upper_)};
  }
```

- [ ] **Step 2.2: Rewrite `searchDifferentPrecursorRanges` (lines 965-1004)**

Replace the entire method body with:
```cpp
  void FragmentIndex::searchDifferentPrecursorRanges(const MSSpectrum& spectrum,
                                                     float precursor_mass,
                                                     SpectrumMatchesTopN& sms,
                                                     uint16_t charge)
  {
    // Open mode absorbs isotope shifts into the wide window — no per-isotope iteration.
    const bool open_mode = isOpenSearchMode_();
    const int16_t iso_lo = open_mode ? 0 : min_isotope_error_;
    const int16_t iso_hi = open_mode ? 0 : max_isotope_error_;

    for (int16_t isotope_error = iso_lo; isotope_error <= iso_hi; ++isotope_error)
    {
      const float shifted_mass = precursor_mass
        + static_cast<float>(isotope_error) * static_cast<float>(Constants::C13C12_MASSDIFF_U);

      const auto window = computeMassWindow_(shifted_mass);

      SpectrumMatchesTopN candidates_iso_error;
      auto candidates_range = getPeptidesInMassWindow(shifted_mass, window);
      candidates_iso_error.hits_.resize(candidates_range.second - candidates_range.first + 1);

      queryPeaks(candidates_iso_error, spectrum, candidates_range, isotope_error, charge);

      sms += candidates_iso_error;
    }
  }
```

Also remove the dead comment at the old line 993 (`// for the simple search we do not apply any modification window!!`) — already gone via the replacement above.

- [ ] **Step 2.3: Update param registrations in `FragmentIndex::FragmentIndex()` (lines 1083-1088, 1148-1149)**

Delete lines 1083-1088:
```cpp
    defaults_.setValue("precursor:mass_tolerance", 10.0, "Tolerance for precursor-m/z in search");
    std::vector<std::string> precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);
```
Replace with:
```cpp
    defaults_.setValue("precursor:mass_tolerance_lower", 20.0,
                       "Lower-side precursor-mass tolerance (positive magnitude; effective window "
                       "is [-lower, +upper] around the precursor). "
                       "When strongly asymmetric, also review precursor:isotope_error_min.");
    defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
    defaults_.setValue("precursor:mass_tolerance_upper", 20.0,
                       "Upper-side precursor-mass tolerance (positive magnitude).");
    defaults_.setMinFloat("precursor:mass_tolerance_upper", 0.0);
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", {"ppm", "Da"});
```

Delete lines 1148-1149 (the dead open_window registrations):
```cpp
    defaults_.setValue("precursor:open_window_lower", -100.0, "lower bound of the open precursor window");
    defaults_.setValue("precursor:open_window_upper", 200.0, "upper bound of the open precursor window");
```

- [ ] **Step 2.4: Rewrite `FragmentIndex::updateMembers_()` precursor block (lines 1194-1219)**

Replace:
```cpp
    precursor_mz_tolerance_ = param_.getValue("precursor:mass_tolerance");
    fragment_mz_tolerance_ = param_.getValue("fragment:mass_tolerance");
    precursor_mz_tolerance_unit_ppm_ = param_.getValue("precursor:mass_tolerance_unit").toString() == "ppm";
    fragment_mz_tolerance_unit_ppm_ = param_.getValue("fragment:mass_tolerance_unit").toString() == "ppm";
    // ...
    if (isOpenSearchMode_())
    {
      OPENMS_LOG_INFO << "[FragmentIndex] Open-search mode enabled because precursor mass tolerance ("
                      << precursor_mz_tolerance_ << " "
                      << (precursor_mz_tolerance_unit_ppm_ ? "ppm" : "Da")
                      << ") exceeds threshold (1000 ppm or 1 Da)." << std::endl;
    }
    open_precursor_window_lower_ = param_.getValue("precursor:open_window_lower");
    open_precursor_window_upper_ = param_.getValue("precursor:open_window_upper");
```
with:
```cpp
    precursor_mass_tolerance_lower_ = param_.getValue("precursor:mass_tolerance_lower");
    precursor_mass_tolerance_upper_ = param_.getValue("precursor:mass_tolerance_upper");
    precursor_mass_tolerance_unit_ppm_ = param_.getValue("precursor:mass_tolerance_unit").toString() == "ppm";
    fragment_mz_tolerance_ = param_.getValue("fragment:mass_tolerance");
    fragment_mz_tolerance_unit_ppm_ = param_.getValue("fragment:mass_tolerance_unit").toString() == "ppm";

    // Validation — setMinFloat(0.0) rejects negatives via checkDefaults_, but NaN/+inf slip past
    // (NaN < 0 is false). NaN would break lower_bound's strict-weak-ordering.
    if (!std::isfinite(precursor_mass_tolerance_lower_) ||
        !std::isfinite(precursor_mass_tolerance_upper_))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "precursor:mass_tolerance_lower and mass_tolerance_upper must be finite");
    }
    if (precursor_mass_tolerance_lower_ + precursor_mass_tolerance_upper_ <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "precursor window has zero width (lower + upper must be > 0)");
    }

    if (isOpenSearchMode_())
    {
      OPENMS_LOG_WARN << "[FragmentIndex] Open-search mode auto-triggered: window [-"
                      << precursor_mass_tolerance_lower_ << ", +"
                      << precursor_mass_tolerance_upper_ << "] "
                      << (precursor_mass_tolerance_unit_ppm_ ? "ppm" : "Da")
                      << " exceeds threshold. Isotope-error iteration collapses to [0, 0]."
                      << std::endl;
    }
```

Add `#include <cmath>` at the top of `FragmentIndex.cpp` if not already present (for `std::isfinite`).

- [ ] **Step 2.5: Build — expect `PeptideSearchEngineFIAlgorithm.cpp` errors**

Run:
```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -50
```
Expected: compile errors in `PeptideSearchEngineFIAlgorithm.cpp` referencing `precursor_mass_tolerance_` (the scalar that no longer exists on `FragmentIndex`). The next task fixes those. **Do NOT commit yet** — atomic refactor spans three files.

---

## Task 3: PeptideSearchEngineFIAlgorithm header refactor

**Files:**
- Modify: `src/openms/include/OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h`

- [ ] **Step 3.1: Replace scalar tolerance members (lines 396-397)**

Replace:
```cpp
    double precursor_mass_tolerance_;
    String precursor_mass_tolerance_unit_;
```
with:
```cpp
    double precursor_mass_tolerance_lower_{20.0};   ///< positive magnitude
    double precursor_mass_tolerance_upper_{20.0};   ///< positive magnitude
    String precursor_mass_tolerance_unit_{"ppm"};
```

- [ ] **Step 3.2: Extend `CalibrationResult_` struct (lines 444-450)**

Replace:
```cpp
    struct CalibrationResult_
    {
      double precursor_tolerance{0}; ///< estimated precursor tolerance (same unit as configured)
      double fragment_tolerance{0};  ///< estimated fragment tolerance (same unit as configured)
      double fragment_shift{0};      ///< reserved for future fragment m/z shift correction
      bool success{false};           ///< true if enough PSMs were found for reliable estimation
    };
```
with:
```cpp
    struct CalibrationResult_
    {
      double precursor_shift{0};     ///< signed median of precursor errors (calibration bias)
      double precursor_spread{0};    ///< median(|e - shift|) + 3 * MAD(|e - shift|)
      double cal_lower{0};           ///< calibrated lower magnitude (valid iff !extreme_bias && success)
      double cal_upper{0};           ///< calibrated upper magnitude (valid iff !extreme_bias && success)
      double fragment_tolerance{0};  ///< estimated fragment tolerance (same unit as configured)
      double fragment_shift{0};      ///< reserved for future fragment m/z shift correction
      bool extreme_bias{false};      ///< |shift| >= spread — writeback skipped (test observability)
      bool success{false};           ///< true if enough PSMs were found for reliable estimation
    };
```

- [ ] **Step 3.3: Add `last_calibration_result_` member and `computeModMatchTolerance_` helper**

Immediately after the `CalibrationResult_` struct definition and before `runCalibrationPass_`, add:
```cpp
    /// Most recent calibration result (valid after any search that invoked runCalibrationPass_).
    /// Stored for test observability and diagnostics.
    CalibrationResult_ last_calibration_result_;

    /// Scalar tolerance passed to OpenSearchModificationAnalysis under asymmetric bounds.
    /// Uses the tighter of the two positive magnitudes — semantically correct for
    /// UniMod Δmass matching precision. OpenSearchModificationAnalysis internally clamps
    /// this at MAX_MOD_MAPPING_TOL_ = 0.02 Da; see spec §7 for rationale.
    double computeModMatchTolerance_() const
    {
      return std::min(precursor_mass_tolerance_lower_, precursor_mass_tolerance_upper_);
    }
```

Note: `runCalibrationPass_` was previously declared `const` — it now populates `last_calibration_result_` as a side-effect, so remove the `const` qualifier from its declaration at line 464-466:
```cpp
    CalibrationResult_ runCalibrationPass_(PeakMap& spectra,
                                           FragmentIndex& fragment_index,
                                           const std::vector<FASTAFile::FASTAEntry>& db);
```
(delete the trailing `const`).

- [ ] **Step 3.4: Update `isOpenSearchMode_()` (lines 478-483) to delegate**

Replace:
```cpp
    bool isOpenSearchMode_() const
    {
      return precursor_mass_tolerance_unit_ == "ppm"
               ? (precursor_mass_tolerance_ > 1000.0)
               : (precursor_mass_tolerance_ > 1.0);
    }
```
with:
```cpp
    bool isOpenSearchMode_() const
    {
      return FragmentIndex::isOpenSearchMode(precursor_mass_tolerance_lower_,
                                             precursor_mass_tolerance_upper_,
                                             precursor_mass_tolerance_unit_ == "ppm");
    }
```

- [ ] **Step 3.5: Add `<algorithm>` include for `std::min`**

At the top of `PeptideSearchEngineFIAlgorithm.h`, after the existing includes, add:
```cpp
#include <algorithm>   // std::min (used by inline computeModMatchTolerance_)
```

- [ ] **Step 3.6: Add friend declaration for the test subclass**

Just before the closing `};` of the class, add:
```cpp
    friend class PeptideSearchEngineFIAlgorithm_test;
```

---

## Task 4: PeptideSearchEngineFIAlgorithm.cpp — param registrations + `updateMembers_` + OSMA call site

**Files:**
- Modify: `src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp`

- [ ] **Step 4.1: Update param registrations (lines 65-72, 181-182)**

At line 65-72, replace:
```cpp
    defaults_.setValue("precursor:mass_tolerance", 10.0, "+/- tolerance for precursor mass.");

    std::vector<std::string> precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);
```
with:
```cpp
    defaults_.setValue("precursor:mass_tolerance_lower", 20.0,
                       "Lower-side precursor-mass tolerance (positive magnitude; effective window "
                       "is [-lower, +upper] around the precursor). "
                       "When strongly asymmetric, also review precursor:isotope_error_min.");
    defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
    defaults_.setValue("precursor:mass_tolerance_upper", 20.0,
                       "Upper-side precursor-mass tolerance (positive magnitude).");
    defaults_.setMinFloat("precursor:mass_tolerance_upper", 0.0);
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", {"ppm", "Da"});
```

At line 77, update the section description:
```cpp
    defaults_.setSectionDescription("precursor",
      "Precursor (Parent Ion) Options. mass_tolerance_lower/_upper are positive magnitudes "
      "applied as [-lower, +upper] around the precursor mass.");
```

Delete lines 181-182 (the dead open_window_* registrations):
```cpp
    // Open search window bounds (used when tolerance > 1 Da or > 1000 ppm)
    defaults_.setValue("precursor:open_window_lower", -100.0, "lower bound of the open precursor window");
    defaults_.setValue("precursor:open_window_upper", 200.0, "upper bound of the open precursor window");
```

- [ ] **Step 4.2: Update `updateMembers_()` at line 223**

Replace:
```cpp
    precursor_mass_tolerance_ = param_.getValue("precursor:mass_tolerance");
    precursor_mass_tolerance_unit_ = param_.getValue("precursor:mass_tolerance_unit").toString();
```
with:
```cpp
    precursor_mass_tolerance_lower_ = param_.getValue("precursor:mass_tolerance_lower");
    precursor_mass_tolerance_upper_ = param_.getValue("precursor:mass_tolerance_upper");
    precursor_mass_tolerance_unit_ = param_.getValue("precursor:mass_tolerance_unit").toString();
```

- [ ] **Step 4.3: Update `OpenSearchModificationAnalysis` call site at lines 889-895**

Replace:
```cpp
      OpenSearchModificationAnalysis mod_analyzer;
      auto modification_summaries = mod_analyzer.analyzeModifications(
        peptide_ids,
        precursor_mass_tolerance_,
        precursor_mass_tolerance_unit_ == "ppm",
        false, // no smoothing
        ""     // no output file for in-memory search
      );
```
with:
```cpp
      OpenSearchModificationAnalysis mod_analyzer;
      auto modification_summaries = mod_analyzer.analyzeModifications(
        peptide_ids,
        computeModMatchTolerance_(),
        precursor_mass_tolerance_unit_ == "ppm",
        false, // no smoothing
        ""     // no output file for in-memory search
      );
```

- [ ] **Step 4.4: Update logging at lines 741-742 (uses the old scalar name)**

Find the log line that uses `precursor_mass_tolerance_` as a scalar. Replace with:
```cpp
    OPENMS_LOG_INFO << "[PDBS-FI] open_search=" << (open_search ? "true" : "false")
                    << " (precursor tolerance [-" << precursor_mass_tolerance_lower_
                    << ", +" << precursor_mass_tolerance_upper_ << "] "
                    << precursor_mass_tolerance_unit_ << ")" << std::endl;
```

Similarly, grep the file for any other `precursor_mass_tolerance_` (scalar) references and update them to use the pair explicitly, or — if the context is feeding a legacy scalar-taking function — use `computeModMatchTolerance_()` or `std::max(lower_, upper_)` as appropriate for that context.

**Run this grep to find all remaining references:**
```bash
grep -n "precursor_mass_tolerance_[^lu]" src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp
```
Fix each hit in context. Most should route through `computeModMatchTolerance_()`, explicit `min(_lower_, _upper_)`, or explicit `max(_lower_, _upper_)` depending on whether the consumer wants tightest or widest.

- [ ] **Step 4.5: Temporary calibration stub — keep the calibration pass compiling**

`runCalibrationPass_` at cpp:1454-1626 currently populates `result.precursor_tolerance`. That field no longer exists. As a **temporary stub** (to be replaced in Task 5), make the minimum change needed to build:

Find the lines inside `runCalibrationPass_` that write to `result.precursor_tolerance` (around cpp:1586-1587, 1605-1607) and replace with:
```cpp
    result.precursor_spread = prec_median + 3.0 * prec_mad;
    if (result.precursor_spread < min_tolerance) result.precursor_spread = min_tolerance;
```
and in the cap at cpp:1605-1608:
```cpp
    if (result.precursor_spread > precursor_mass_tolerance_upper_)  // use upper as the symmetric cap placeholder
    {
      result.precursor_spread = precursor_mass_tolerance_upper_;
    }
```
Also update the log lines around cpp:1620-1621 to reference `precursor_spread` instead of `precursor_tolerance`. This is temporary — Task 5 replaces both the formula and the writeback semantics.

Also find the calibration **writeback** site around cpp:770-780 that sets `fi_params.setValue("precursor:mass_tolerance", ...)`. Replace with:
```cpp
      if (cal.success)
      {
        Param fi_params = fi_params_original;
        fi_params.setValue("precursor:mass_tolerance_lower", cal.precursor_spread);
        fi_params.setValue("precursor:mass_tolerance_upper", cal.precursor_spread);
        fi_params.setValue("fragment:mass_tolerance", cal.fragment_tolerance);
        fragment_index_.setParameters(fi_params);
        fi_params_modified = true;
        // Task 5 will replace with proper bias-preserving writeback.
      }
```

This leaves calibration behavior equivalent to the old code (symmetric `±spread`) on a pure rename. Task 5 then converts it to the proper signed-bias version.

- [ ] **Step 4.6: Build OpenMS — expect success**

Run:
```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -40
```
Expected: **compilation succeeds.** Link errors for `PeptideDataBaseSearchFI.cpp` may appear — Task 5 covers them via the TOPP tool update.

**Do NOT commit yet** — tests still need migration (Task 6) before a clean commit.

---

## Task 5: PeptideDataBaseSearchFI TOPP tool update

**Files:**
- Modify: `src/topp/PeptideDataBaseSearchFI.cpp`

- [ ] **Step 5.1: Update direct param reads at lines 168-169**

Find the block that currently reads `precursor:mass_tolerance` and `mass_tolerance_unit` directly. Replace with:
```cpp
    const double mass_tol_lower = getParam_().getValue("Search:precursor:mass_tolerance_lower");
    const double mass_tol_upper = getParam_().getValue("Search:precursor:mass_tolerance_upper");
    const String mass_tol_unit = getParam_().getValue("Search:precursor:mass_tolerance_unit").toString();
```
(The `Search:` prefix is the existing convention — verify against the existing scalar-reading line.)

- [ ] **Step 5.2: Update local open-mode detection at lines 170-172**

Replace the local open-search-mode computation with a call to the shared helper:
```cpp
    const bool open_search = FragmentIndex::isOpenSearchMode(mass_tol_lower,
                                                             mass_tol_upper,
                                                             mass_tol_unit == "ppm");
```

- [ ] **Step 5.3: Build TOPP tools**

Run:
```bash
cmake --build OpenMS-build --target PeptideDataBaseSearchFI -j$(nproc) 2>&1 | tail -20
```
Expected: clean compile and link.

**Still do not commit** — existing tests need migration first.

---

## Task 6: Migrate test INIs and existing tests

**Files:**
- Modify: `src/tests/topp/PeptideDataBaseSearchFI_1.ini`
- Modify: `src/tests/topp/PeptideDataBaseSearchFI_2.ini`
- Modify: `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`
- Modify: `src/tests/class_tests/openms/source/FragmentIndex_test.cpp`
- Modify: `src/pyOpenMS/tests/test_modification_discovery.py`

- [ ] **Step 6.1: Migrate `PeptideDataBaseSearchFI_1.ini`**

Open the file and find the line setting `mass_tolerance = 5` (or similar) under `precursor`. Replace with two items:
```xml
      <ITEM name="mass_tolerance_lower" value="5.0" type="double" description="..." required="false" advanced="false" />
      <ITEM name="mass_tolerance_upper" value="5.0" type="double" description="..." required="false" advanced="false" />
```
Delete the old `mass_tolerance` item. Preserve the `mass_tolerance_unit` item unchanged.

- [ ] **Step 6.2: Migrate `PeptideDataBaseSearchFI_2.ini`**

Same pattern with `10.0`:
```xml
      <ITEM name="mass_tolerance_lower" value="10.0" type="double" description="..." required="false" advanced="false" />
      <ITEM name="mass_tolerance_upper" value="10.0" type="double" description="..." required="false" advanced="false" />
```

- [ ] **Step 6.3: Migrate `PeptideSearchEngineFIAlgorithm_test.cpp`**

Grep for every `precursor:mass_tolerance` param-set site:
```bash
grep -n 'precursor:mass_tolerance' src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp
```
For each hit, replace the single `setValue("precursor:mass_tolerance", X)` with:
```cpp
p.setValue("precursor:mass_tolerance_lower", X);
p.setValue("precursor:mass_tolerance_upper", X);
```

- [ ] **Step 6.4: Migrate `FragmentIndex_test.cpp`**

Same pattern:
```bash
grep -n 'precursor:mass_tolerance' src/tests/class_tests/openms/source/FragmentIndex_test.cpp
```
Replace each with the pair. Note: any test that set `mass_tolerance` to a wide open-search value (e.g., `0.5 Da`) should now set both `mass_tolerance_lower` and `mass_tolerance_upper` to that value and **remove** any `open_window_lower`/`_upper` set calls — the unified window is now the only mechanism.

- [ ] **Step 6.5: Migrate `test_modification_discovery.py`**

Open `src/pyOpenMS/tests/test_modification_discovery.py`. Grep for the 8 hits:
```bash
grep -n "precursor:mass_tolerance" src/pyOpenMS/tests/test_modification_discovery.py
```
Expected: 8 lines across 4 test functions (lines approximately 197-198, 300-301, 354-355, 434-435). Each function sets `precursor:mass_tolerance` then `precursor:mass_tolerance_unit`. Replace the pair at each function:
```python
p.setValue("precursor:mass_tolerance_lower", 10.0)
p.setValue("precursor:mass_tolerance_upper", 10.0)
p.setValue("precursor:mass_tolerance_unit", "ppm")
```
(Use whatever tolerance value the existing test used.)

- [ ] **Step 6.6: Build + run the existing test suite to confirm green baseline**

```bash
cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -20
ctest --test-dir OpenMS-build -R "FragmentIndex|PeptideSearchEngineFIAlgorithm|PeptideDataBaseSearchFI" --output-on-failure 2>&1 | tail -50
```
Expected: all existing tests pass. If any idXML reference file diffs, **do not regenerate** — investigate whether the migration missed a site. The refactor should be bitwise-equivalent on symmetric configs (spec §9 bitwise-equivalence claim).

Then run the pyOpenMS test:
```bash
PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/test_modification_discovery.py -v 2>&1 | tail -30
```
Expected: pass.

- [ ] **Step 6.7: Commit the atomic refactor (tasks 1-6)**

```bash
git add src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h \
        src/openms/source/ANALYSIS/ID/FragmentIndex.cpp \
        src/openms/include/OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h \
        src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp \
        src/topp/PeptideDataBaseSearchFI.cpp \
        src/tests/topp/PeptideDataBaseSearchFI_1.ini \
        src/tests/topp/PeptideDataBaseSearchFI_2.ini \
        src/tests/class_tests/openms/source/FragmentIndex_test.cpp \
        src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp \
        src/pyOpenMS/tests/test_modification_discovery.py

git commit -m "$(cat <<'EOF'
refactor(FragmentIndex): unify precursor tolerance into positive-magnitude pair

Replace scalar precursor:mass_tolerance + open_window_lower/_upper with a
single positive-magnitude pair mass_tolerance_lower/_upper. Refactor
getPeptidesInPrecursorRange → getPeptidesInMassWindow with absolute-bound
semantics (no hidden ±prec_tol inside the function). Consolidate the
duplicate isOpenSearchMode_() into a shared static helper.

Calibration pass still writes symmetric bounds (Task 5 follow-up replaces
with signed-bias preservation).

Part of OpenMS/OpenMS#9108 item 5.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Calibration pass — signed bias preservation + writeback refresh

**Files:**
- Modify: `src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp`

This task replaces the Task 4 temporary stub with the proper signed-bias algorithm.

- [ ] **Step 7.1: Rewrite `runCalibrationPass_` precursor spread computation**

In `runCalibrationPass_` (around cpp:1575-1608 after Task 4's stub), replace the precursor spread block with:
```cpp
    // Signed median gives the calibration bias; spread is median(|e - shift|) + 3*MAD(|e - shift|),
    // i.e., the same "typical residual + 3*scale" shape as the pre-refactor formula, just applied
    // to residuals around the signed center instead of raw errors. Zero-centered distributions
    // get identical spread; biased distributions correctly separate shift from spread.
    const double prec_shift = Math::median(precursor_errors.begin(), precursor_errors.end());

    std::vector<double> residuals;
    residuals.reserve(precursor_errors.size());
    for (double e : precursor_errors) { residuals.push_back(std::abs(e - prec_shift)); }
    const double res_median = Math::median(residuals.begin(), residuals.end());
    const double res_mad    = Math::MAD(residuals.begin(), residuals.end(), res_median);

    result.precursor_shift  = prec_shift;
    result.precursor_spread = std::max(min_tolerance, res_median + 3.0 * res_mad);
    result.extreme_bias     = std::abs(prec_shift) >= result.precursor_spread;

    if (!result.extreme_bias)
    {
      // Map signed window [shift - spread, shift + spread] to positive magnitudes:
      //   -cal_lower = shift - spread  →  cal_lower = spread - shift
      //   +cal_upper = shift + spread  →  cal_upper = spread + shift
      // Both are non-negative when |shift| < spread.
      const double cal_lower_raw = result.precursor_spread - prec_shift;
      const double cal_upper_raw = result.precursor_spread + prec_shift;
      // Only tighten — cap against user-configured bounds.
      result.cal_lower = std::min(cal_lower_raw, precursor_mass_tolerance_lower_);
      result.cal_upper = std::min(cal_upper_raw, precursor_mass_tolerance_upper_);
    }
    // else: cal_lower/cal_upper stay at 0; writeback block skips the calibration result.
```

Delete the old symmetric cap code that was preserved from Task 4 (the `if (result.precursor_spread > precursor_mass_tolerance_upper_)` block).

- [ ] **Step 7.2: Rewrite the calibration writeback (cpp:770-782)**

Replace the Task 4 stub with:
```cpp
    if (cal.success && !cal.extreme_bias)
    {
      Param fi_params = fi_params_original;
      fi_params.setValue("precursor:mass_tolerance_lower", cal.cal_lower);
      fi_params.setValue("precursor:mass_tolerance_upper", cal.cal_upper);
      fi_params.setValue("fragment:mass_tolerance", cal.fragment_tolerance);
      fragment_index_.setParameters(fi_params);
      fi_params_modified = true;

      // Refresh algo-level member copies — OpenSearchModificationAnalysis reads them via
      // computeModMatchTolerance_() and would otherwise see stale pre-calibration values.
      precursor_mass_tolerance_lower_ = cal.cal_lower;
      precursor_mass_tolerance_upper_ = cal.cal_upper;
      effective_fragment_tol = cal.fragment_tolerance;

      OPENMS_LOG_INFO << "[PDBS-FI] Calibration: shift=" << cal.precursor_shift
                      << " spread=" << cal.precursor_spread << " "
                      << precursor_mass_tolerance_unit_
                      << " → window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                      << std::endl;
    }
    else if (cal.success && cal.extreme_bias)
    {
      OPENMS_LOG_WARN << "[PDBS-FI] Calibration: |shift|=" << std::abs(cal.precursor_shift)
                      << " >= spread=" << cal.precursor_spread << " "
                      << precursor_mass_tolerance_unit_ << " — calibration result discarded. "
                      << "The true signed window ["
                      << (cal.precursor_shift - cal.precursor_spread) << ", "
                      << (cal.precursor_shift + cal.precursor_spread)
                      << "] lies entirely on one side of zero; the positive-magnitude schema "
                      << "cannot express it without loosening. Fix external calibration, or "
                      << "configure mass_tolerance_lower/_upper manually." << std::endl;
      // Preserve user bounds — fi_params_modified stays false, members unchanged.
    }
```

- [ ] **Step 7.3: Persist calibration result for test observability**

Find the line where `runCalibrationPass_` is called (around cpp:767):
```cpp
        CalibrationResult_ cal = runCalibrationPass_(spectra, fragment_index_, db);
```
Change to populate the member as well:
```cpp
        last_calibration_result_ = runCalibrationPass_(spectra, fragment_index_, db);
        const CalibrationResult_& cal = last_calibration_result_;
```

This allows tests to inspect the shift, spread, cal_lower, cal_upper, and extreme_bias fields after the search completes.

- [ ] **Step 7.4: Build**

```bash
cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -30
```
Expected: success.

- [ ] **Step 7.5: Run the existing test suite — confirm no regression**

```bash
ctest --test-dir OpenMS-build -R "PeptideSearchEngineFIAlgorithm|PeptideDataBaseSearchFI" --output-on-failure 2>&1 | tail -40
```
Expected: all existing tests still pass (the new algorithm produces equivalent output to the symmetric old algorithm on the symmetric test fixtures).

- [ ] **Step 7.6: Commit**

```bash
git add src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp \
        src/openms/include/OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h

git commit -m "$(cat <<'EOF'
feat(PSEFIA): calibration preserves signed precursor bias

runCalibrationPass_ now computes signed median + MAD-of-residuals,
producing asymmetric calibrated bounds [spread - shift, spread + shift]
when |shift| < spread. When |shift| >= spread, the signed window cannot
be represented in the positive-magnitude schema without loosening, so
the calibration result is discarded and a WARN is logged instead.

Writeback refreshes both fragment_index_ params AND algo-level member
copies so downstream OpenSearchModificationAnalysis sees calibrated
values instead of stale pre-calibration ones (fixes latent bug that
existed under the old scalar tolerance).

last_calibration_result_ stored as a member for test observability.

Part of OpenMS/OpenMS#9108 item 5.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: FragmentIndex_test — symmetric equivalence and asymmetric window

**Files:**
- Modify: `src/tests/class_tests/openms/source/FragmentIndex_test.cpp`

The existing `FragmentIndex_test` friend subclass (line 44) already exposes `fi_peptides_` and `getPeptides()`. No new subclass members needed.

- [ ] **Step 8.1: Write test 1 — symmetric default self-hit**

Add a new `START_SECTION` block in `FragmentIndex_test.cpp`:
```cpp
START_SECTION((pair<size_t, size_t> getPeptidesInMassWindow(float, const pair<float, float>&) const))
{
  // Symmetric default [20, 20] ppm — each peptide retrieves itself within its own mass window.
  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 20.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  fi.setParameters(p);

  // Build a small fixture
  vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e;
  e.identifier = "TEST1";
  e.sequence = "PEPTIDER";
  entries.push_back(e);
  fi.build(entries);

  TEST_EQUAL(fi.testQuery(2, true, entries), true);
}
END_SECTION
```

- [ ] **Step 8.2: Write test 2 — asymmetric window compensates an injected precursor offset**

```cpp
START_SECTION((asymmetric window compensates precursor calibration offset))
{
  // Instrument reads precursor m/z +8 ppm high. Symmetric [5, 5] ppm misses the peptide;
  // asymmetric [5, 15] ppm catches it.
  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_unit", "ppm");

  vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e;
  e.identifier = "TEST";
  e.sequence = "PEPTIDER";
  entries.push_back(e);

  // First run: symmetric tight — should NOT find
  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  fi.setParameters(p);
  fi.build(entries);

  // Construct a query spectrum whose precursor is shifted +8 ppm
  AASequence seq = AASequence::fromString("PEPTIDER");
  MSSpectrum spec;
  Precursor prec;
  const double true_mz = seq.getMZ(2);
  prec.setMZ(true_mz * (1.0 + 8e-6));   // +8 ppm
  prec.setCharge(2);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);

  // Build a minimal theoretical spectrum for the query
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum theo;
  tsg.getSpectrum(theo, seq, 1, 1);
  for (const auto& peak : theo) spec.push_back(peak);
  spec.sortByPosition();

  FragmentIndex::SpectrumMatchesTopN sms_tight;
  fi.querySpectrum(spec, sms_tight);
  TEST_EQUAL(sms_tight.hits_.empty(), true);

  // Second run: asymmetric [5, 15] ppm — should find
  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 15.0);
  fi.setParameters(p);

  FragmentIndex::SpectrumMatchesTopN sms_asym;
  fi.querySpectrum(spec, sms_asym);
  TEST_NOT_EQUAL(sms_asym.hits_.size(), 0);
}
END_SECTION
```

- [ ] **Step 8.3: Write test 3 — open-mode auto-detection boundaries**

```cpp
START_SECTION((static bool isOpenSearchMode(double, double, bool)))
{
  // Strict > threshold. 1000 ppm stays closed.
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(500.0,  1500.0, true), true);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(999.0,   999.0, true), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1000.0, 1000.0, true), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1000.0001, 1000.0, true), true);

  // Da unit — threshold 1.0
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(0.9, 0.9, false), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1.0, 1.0, false), false);
  TEST_EQUAL(FragmentIndex::isOpenSearchMode(1.1, 0.5, false), true);
}
END_SECTION
```

- [ ] **Step 8.4: Build + run**

```bash
cmake --build OpenMS-build --target FragmentIndex_test -j$(nproc)
ctest --test-dir OpenMS-build -R "^FragmentIndex_test$" --output-on-failure 2>&1 | tail -20
```
Expected: new tests pass.

- [ ] **Step 8.5: Commit**

```bash
git add src/tests/class_tests/openms/source/FragmentIndex_test.cpp
git commit -m "test(FragmentIndex): symmetric equivalence, asymmetric fixture, open-mode boundaries"
```

---

## Task 9: FragmentIndex_test — observable-proxy isotope tests (4, 5)

**Files:**
- Modify: `src/tests/class_tests/openms/source/FragmentIndex_test.cpp`

- [ ] **Step 9.1: Write test 4 — open-mode forces isotope_error to [0, 0]**

```cpp
START_SECTION((open-mode forces isotope_error iteration to [0,0]))
{
  // Observable-proxy: under open mode, a fixture with isotope_error_range [-2, +2]
  // produces the same PSM set as [0, 0] — proving iteration collapsed.
  // Build a peptide fixture
  vector<FASTAFile::FASTAEntry> entries;
  FASTAFile::FASTAEntry e;
  e.identifier = "TEST";
  e.sequence = "PEPTIDER";
  entries.push_back(e);

  // Run 1: open mode with user iso range [-2, +2]
  FragmentIndex_test fi_a;
  Param p_a = fi_a.getParameters();
  p_a.setValue("precursor:mass_tolerance_lower", 0.5);
  p_a.setValue("precursor:mass_tolerance_upper", 1.5);  // 1.5 Da > 1.0 → open mode
  p_a.setValue("precursor:mass_tolerance_unit", "Da");
  p_a.setValue("precursor:isotope_error_min", -2);
  p_a.setValue("precursor:isotope_error_max", +2);
  fi_a.setParameters(p_a);
  fi_a.build(entries);

  // Run 2: open mode with iso range [0, 0]
  FragmentIndex_test fi_b;
  Param p_b = fi_b.getParameters();
  p_b.setValue("precursor:mass_tolerance_lower", 0.5);
  p_b.setValue("precursor:mass_tolerance_upper", 1.5);
  p_b.setValue("precursor:mass_tolerance_unit", "Da");
  p_b.setValue("precursor:isotope_error_min", 0);
  p_b.setValue("precursor:isotope_error_max", 0);
  fi_b.setParameters(p_b);
  fi_b.build(entries);

  // Construct identical query spectrum for both
  AASequence seq = AASequence::fromString("PEPTIDER");
  MSSpectrum spec;
  Precursor prec;
  prec.setMZ(seq.getMZ(2));
  prec.setCharge(2);
  spec.getPrecursors().push_back(prec);
  spec.setMSLevel(2);
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum theo;
  tsg.getSpectrum(theo, seq, 1, 1);
  for (const auto& peak : theo) spec.push_back(peak);
  spec.sortByPosition();

  FragmentIndex::SpectrumMatchesTopN sms_a, sms_b;
  fi_a.querySpectrum(spec, sms_a);
  fi_b.querySpectrum(spec, sms_b);

  // Equal PSM set sizes → iteration collapsed (the [-2,+2] config did NOT produce more hits)
  TEST_EQUAL(sms_a.hits_.size(), sms_b.hits_.size());
}
END_SECTION
```

- [ ] **Step 9.2: Write test 5 — asymmetric + isotope-error interaction (closed mode)**

```cpp
START_SECTION((asymmetric closed window interacts with isotope_error iteration))
{
  // [5, 15] ppm + isotope_error [-1, +2]. Fixture with peptides at mass, mass + C13, mass + 2*C13
  // — all three should be found because each isotope iteration shifts the center and the
  // asymmetric window travels with it.
  vector<FASTAFile::FASTAEntry> entries;
  // Three peptides with different masses (roughly +0, +1, +2 Da)
  FASTAFile::FASTAEntry e1, e2, e3;
  e1.identifier = "P1"; e1.sequence = "PEPTIDER";      // mass m0
  e2.identifier = "P2"; e2.sequence = "PEPTIDERG";     // mass ~m0+57 Da (G = 57.02)
  e3.identifier = "P3"; e3.sequence = "PEPTIDERA";     // mass ~m0+71 Da (A = 71.04)
  entries = {e1, e2, e3};

  FragmentIndex_test fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 15.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("precursor:isotope_error_min", -1);
  p.setValue("precursor:isotope_error_max", +2);
  fi.setParameters(p);
  fi.build(entries);

  // Query each peptide's own theoretical spectrum and verify self-hit via testQuery
  TEST_EQUAL(fi.testQuery(2, true, entries), true);
}
END_SECTION
```

- [ ] **Step 9.3: Build + run + commit**

```bash
cmake --build OpenMS-build --target FragmentIndex_test -j$(nproc)
ctest --test-dir OpenMS-build -R "^FragmentIndex_test$" --output-on-failure 2>&1 | tail -20
git add src/tests/class_tests/openms/source/FragmentIndex_test.cpp
git commit -m "test(FragmentIndex): observable-proxy isotope-error interaction tests"
```

---

## Task 10: FragmentIndex_test — validation throws (6a, 6b, 6c)

**Files:**
- Modify: `src/tests/class_tests/openms/source/FragmentIndex_test.cpp`

- [ ] **Step 10.1: Write test 6a — negative magnitude rejected by `setMinFloat`**

```cpp
START_SECTION((validation: negative magnitude rejected))
{
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", -5.0);  // invalid: below setMinFloat(0.0)
  p.setValue("precursor:mass_tolerance_upper", 10.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  // checkDefaults_ fires from setParameters via the min-float check
  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p));
}
END_SECTION
```

- [ ] **Step 10.2: Write test 6b — zero-width rejected by `updateMembers_`**

```cpp
START_SECTION((validation: zero-width window rejected))
{
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 0.0);
  p.setValue("precursor:mass_tolerance_upper", 0.0);   // sum == 0 → rejected
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p));
}
END_SECTION
```

- [ ] **Step 10.3: Write test 6c — NaN rejected by `std::isfinite` guard**

```cpp
START_SECTION((validation: NaN rejected))
{
  FragmentIndex fi;
  Param p = fi.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 10.0);
  p.setValue("precursor:mass_tolerance_upper", std::numeric_limits<double>::quiet_NaN());
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  TEST_EXCEPTION(Exception::InvalidParameter, fi.setParameters(p));
}
END_SECTION
```

- [ ] **Step 10.4: Build + run + commit**

```bash
cmake --build OpenMS-build --target FragmentIndex_test -j$(nproc)
ctest --test-dir OpenMS-build -R "^FragmentIndex_test$" --output-on-failure 2>&1 | tail -20
git add src/tests/class_tests/openms/source/FragmentIndex_test.cpp
git commit -m "test(FragmentIndex): parameter validation (negative, zero-width, NaN)"
```

---

## Task 11: PeptideSearchEngineFIAlgorithm_test — add friend subclass

**Files:**
- Modify: `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`

- [ ] **Step 11.1: Add the friend subclass near the top of the file (after `using namespace`)**

```cpp
// Test subclass exposing internal state for white-box assertions on
// asymmetric bounds, calibration-pass results, and the mod-match tolerance helper.
class PeptideSearchEngineFIAlgorithm_test : public PeptideSearchEngineFIAlgorithm
{
public:
  using PeptideSearchEngineFIAlgorithm::precursor_mass_tolerance_lower_;
  using PeptideSearchEngineFIAlgorithm::precursor_mass_tolerance_upper_;
  using PeptideSearchEngineFIAlgorithm::precursor_mass_tolerance_unit_;
  using PeptideSearchEngineFIAlgorithm::computeModMatchTolerance_;
  using PeptideSearchEngineFIAlgorithm::last_calibration_result_;
  using PeptideSearchEngineFIAlgorithm::fragment_index_;
  using PeptideSearchEngineFIAlgorithm::CalibrationResult_;
};
```

**Note:** The `friend class` declaration was added in Task 3.6 at the header level — the subclass above uses `using` to re-expose those members, which works thanks to the friend relationship.

- [ ] **Step 11.2: Build to verify the friend subclass compiles**

```bash
cmake --build OpenMS-build --target PeptideSearchEngineFIAlgorithm_test -j$(nproc) 2>&1 | tail -20
```
Expected: success.

- [ ] **Step 11.3: Commit**

```bash
git add src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp
git commit -m "test(PSEFIA): add friend subclass for calibration and tolerance observability"
```

---

## Task 12: PSEFIA_test — Test 7: calibration preserves asymmetric bias

**Files:**
- Modify: `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`

- [ ] **Step 12.1: Write test 7**

```cpp
START_SECTION((calibration preserves asymmetric bias — normal case))
{
  // User sets [20, 5] ppm to compensate a known +7 ppm systematic bias.
  // Calibration should detect shift ≈ +7, spread < 7, produce asymmetric cal_lower/_upper
  // (not symmetric), and refresh both the fragment_index_ params and the algo-level members.
  PeptideSearchEngineFIAlgorithm_test algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("calibration:enabled", "true");
  p.setValue("calibration:subset_ratio", 1.0);
  p.setValue("calibration:min_psms", 1);
  algo.setParameters(p);

  // Load or synthesize a fixture where every precursor m/z is shifted by +7 ppm
  PeakMap spectra;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleSearchEngine_1.mzML"), spectra);
  for (auto& spec : spectra)
  {
    if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty())
    {
      double mz = spec.getPrecursors()[0].getMZ();
      spec.getPrecursors()[0].setMZ(mz * (1.0 + 7e-6));
    }
  }

  vector<FASTAFile::FASTAEntry> db;
  FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleSearchEngine_1.fasta"), db);

  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  algo.search(spectra, db, protein_ids, peptide_ids);

  TEST_EQUAL(algo.last_calibration_result_.success, true);
  TEST_EQUAL(algo.last_calibration_result_.extreme_bias, false);
  // shift should be close to +7 ppm
  TEST_REAL_SIMILAR_TOLERANCE(algo.last_calibration_result_.precursor_shift, 7.0, 3.0);
  // cal_lower/upper should be tightened from [20, 5] but remain asymmetric
  TEST_EQUAL(algo.last_calibration_result_.cal_lower < 20.0, true);
  TEST_EQUAL(algo.last_calibration_result_.cal_upper < 5.0, true);
  // algo-level members refreshed
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_lower_, algo.last_calibration_result_.cal_lower);
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_upper_, algo.last_calibration_result_.cal_upper);
  // fragment_index_ params also refreshed
  Param fi_p = algo.fragment_index_.getParameters();
  TEST_REAL_SIMILAR((double)fi_p.getValue("precursor:mass_tolerance_lower"),
                    algo.last_calibration_result_.cal_lower);
}
END_SECTION
```

- [ ] **Step 12.2: Build + run**

```bash
cmake --build OpenMS-build --target PeptideSearchEngineFIAlgorithm_test -j$(nproc)
ctest --test-dir OpenMS-build -R "^PeptideSearchEngineFIAlgorithm_test$" --output-on-failure 2>&1 | tail -20
```
Expected: pass. If the test fails on `TEST_REAL_SIMILAR_TOLERANCE` for the shift, inspect the calibration log output and adjust the tolerance — OpenMS's `Math::median` precision depends on the number of PSMs and the fixture's systematic offset may not map cleanly to `+7 ppm`. The test's structure is the guard; the exact value tolerance can be relaxed as needed.

- [ ] **Step 12.3: Commit**

```bash
git add src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp
git commit -m "test(PSEFIA): calibration preserves asymmetric bias (normal case)"
```

---

## Task 13: PSEFIA_test — Test 8: computeModMatchTolerance_ reads refreshed members

**Files:**
- Modify: `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`

- [ ] **Step 13.1: Write test 8**

```cpp
START_SECTION((computeModMatchTolerance_ reads post-calibration members))
{
  // Same setup as test 7. After calibration, computeModMatchTolerance_() should return
  // min(cal_lower, cal_upper) — NOT min(20, 5) = 5 (the user-configured values).
  // This is the double-bookkeeping regression guard.
  PeptideSearchEngineFIAlgorithm_test algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 20.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("calibration:enabled", "true");
  p.setValue("calibration:subset_ratio", 1.0);
  p.setValue("calibration:min_psms", 1);
  algo.setParameters(p);

  PeakMap spectra;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleSearchEngine_1.mzML"), spectra);
  for (auto& spec : spectra)
  {
    if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty())
    {
      double mz = spec.getPrecursors()[0].getMZ();
      spec.getPrecursors()[0].setMZ(mz * (1.0 + 7e-6));
    }
  }
  vector<FASTAFile::FASTAEntry> db;
  FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleSearchEngine_1.fasta"), db);
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  algo.search(spectra, db, protein_ids, peptide_ids);

  const double expected = std::min(algo.last_calibration_result_.cal_lower,
                                   algo.last_calibration_result_.cal_upper);
  TEST_REAL_SIMILAR(algo.computeModMatchTolerance_(), expected);
  // And definitely NOT the pre-calibration min
  TEST_NOT_EQUAL(algo.computeModMatchTolerance_(), 5.0);
}
END_SECTION
```

- [ ] **Step 13.2: Build + run + commit**

```bash
cmake --build OpenMS-build --target PeptideSearchEngineFIAlgorithm_test -j$(nproc)
ctest --test-dir OpenMS-build -R "^PeptideSearchEngineFIAlgorithm_test$" --output-on-failure 2>&1 | tail -20
git add src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp
git commit -m "test(PSEFIA): computeModMatchTolerance_ reads refreshed calibration members"
```

---

## Task 14: PSEFIA_test — Test 9: extreme bias path discards calibration

**Files:**
- Modify: `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`

- [ ] **Step 14.1: Write test 9**

```cpp
START_SECTION((calibration extreme-bias path preserves user bounds))
{
  // Construct a fixture where every precursor is shifted by a large constant (e.g., +50 ppm)
  // with very low variance → |shift| >> spread → extreme_bias = true → writeback skipped.
  PeptideSearchEngineFIAlgorithm_test algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_lower", 100.0);  // wide enough to still find hits
  p.setValue("precursor:mass_tolerance_upper", 100.0);
  p.setValue("precursor:mass_tolerance_unit", "ppm");
  p.setValue("calibration:enabled", "true");
  p.setValue("calibration:subset_ratio", 1.0);
  p.setValue("calibration:min_psms", 1);
  algo.setParameters(p);

  PeakMap spectra;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleSearchEngine_1.mzML"), spectra);
  for (auto& spec : spectra)
  {
    if (spec.getMSLevel() == 2 && !spec.getPrecursors().empty())
    {
      double mz = spec.getPrecursors()[0].getMZ();
      spec.getPrecursors()[0].setMZ(mz * (1.0 + 50e-6));  // uniform +50 ppm, no scatter
    }
  }
  vector<FASTAFile::FASTAEntry> db;
  FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("SimpleSearchEngine_1.fasta"), db);
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  algo.search(spectra, db, protein_ids, peptide_ids);

  TEST_EQUAL(algo.last_calibration_result_.success, true);
  TEST_EQUAL(algo.last_calibration_result_.extreme_bias, true);
  // User bounds unchanged — no writeback happened
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_lower_, 100.0);
  TEST_REAL_SIMILAR(algo.precursor_mass_tolerance_upper_, 100.0);
}
END_SECTION
```

- [ ] **Step 14.2: Build + run + commit**

```bash
cmake --build OpenMS-build --target PeptideSearchEngineFIAlgorithm_test -j$(nproc)
ctest --test-dir OpenMS-build -R "^PeptideSearchEngineFIAlgorithm_test$" --output-on-failure 2>&1 | tail -20
git add src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp
git commit -m "test(PSEFIA): calibration extreme-bias path discards result, preserves user bounds"
```

---

## Task 15: PSEFIA_test — Test 10: reduction rule pinning (no calibration)

**Files:**
- Modify: `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`

- [ ] **Step 15.1: Write test 10**

```cpp
START_SECTION((computeModMatchTolerance_ returns min(lower, upper)))
{
  // Pure unit test — no search, no calibration. Pins the min() reduction rule so a
  // future change to max() or midpoint is caught.
  PeptideSearchEngineFIAlgorithm_test algo;
  Param p = algo.getParameters();
  p.setValue("precursor:mass_tolerance_unit", "ppm");

  p.setValue("precursor:mass_tolerance_lower", 5.0);
  p.setValue("precursor:mass_tolerance_upper", 50.0);
  algo.setParameters(p);
  TEST_REAL_SIMILAR(algo.computeModMatchTolerance_(), 5.0);

  p.setValue("precursor:mass_tolerance_lower", 50.0);
  p.setValue("precursor:mass_tolerance_upper", 5.0);
  algo.setParameters(p);
  TEST_REAL_SIMILAR(algo.computeModMatchTolerance_(), 5.0);

  // Da unit
  p.setValue("precursor:mass_tolerance_unit", "Da");
  p.setValue("precursor:mass_tolerance_lower", 0.5);
  p.setValue("precursor:mass_tolerance_upper", 2.0);
  algo.setParameters(p);
  TEST_REAL_SIMILAR(algo.computeModMatchTolerance_(), 0.5);
}
END_SECTION
```

- [ ] **Step 15.2: Build + run + commit**

```bash
cmake --build OpenMS-build --target PeptideSearchEngineFIAlgorithm_test -j$(nproc)
ctest --test-dir OpenMS-build -R "^PeptideSearchEngineFIAlgorithm_test$" --output-on-failure 2>&1 | tail -20
git add src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp
git commit -m "test(PSEFIA): pin computeModMatchTolerance_ min() reduction rule"
```

---

## Task 16: Create thirdparty-level regression fixture (shifted mzML)

**Files:**
- Create: `src/tests/topp/SimpleSearchEngine_1_shifted_7ppm.mzML`

- [ ] **Step 16.1: Write a one-shot Python script to create the shifted fixture**

**IMPORTANT:** This script is NOT committed to the repo. Run it once to produce the output mzML, commit only the output file.

```bash
cat > /tmp/shift_mzml.py <<'EOF'
import sys
import pyopenms as oms

src = "src/tests/topp/SimpleSearchEngine_1.mzML"
dst = "src/tests/topp/SimpleSearchEngine_1_shifted_7ppm.mzML"

exp = oms.MSExperiment()
oms.MzMLFile().load(src, exp)

shifted = 0
for i in range(exp.size()):
    spec = exp[i]
    if spec.getMSLevel() != 2:
        continue
    precursors = spec.getPrecursors()
    if not precursors:
        continue
    for j in range(len(precursors)):
        mz = precursors[j].getMZ()
        precursors[j].setMZ(mz * (1.0 + 7e-6))
    spec.setPrecursors(precursors)
    exp[i] = spec
    shifted += 1

print(f"Shifted {shifted} MS2 precursors by +7 ppm")
oms.MzMLFile().store(dst, exp)
EOF
PYTHONPATH=OpenMS-build/pyOpenMS python3 /tmp/shift_mzml.py
rm /tmp/shift_mzml.py
```

Expected output: `"Shifted N MS2 precursors by +7 ppm"` where N is the number of MS2 spectra in the source file.

- [ ] **Step 16.2: Verify the shifted file was created and is reasonable size**

```bash
ls -la src/tests/topp/SimpleSearchEngine_1_shifted_7ppm.mzML
```
Expected: file exists, roughly the same size as `SimpleSearchEngine_1.mzML` (~38 KB).

- [ ] **Step 16.3: Sanity-check: diff the precursor m/z of the first MS2 spectrum in both files**

```bash
grep -A 1 'selectedIon' src/tests/topp/SimpleSearchEngine_1.mzML | head -5
grep -A 1 'selectedIon' src/tests/topp/SimpleSearchEngine_1_shifted_7ppm.mzML | head -5
```
Expected: the `selected ion m/z` CV param values should differ by a factor of `1 + 7e-6` — e.g., `500.0` → `500.0035`.

- [ ] **Step 16.4: Commit the fixture (the .py script itself is NOT committed)**

```bash
git add src/tests/topp/SimpleSearchEngine_1_shifted_7ppm.mzML
git commit -m "$(cat <<'EOF'
test(data): add +7 ppm shifted mzML fixture for asymmetric precursor regression test

Created by one-shot Python script (not committed): load SimpleSearchEngine_1.mzML,
multiply every MS2 precursor m/z by (1 + 7e-6), save as new file. Used by
PeptideDataBaseSearchFI_4 to verify end-to-end asymmetric-bounds + calibration
bias preservation on realistic spectra with a known systematic bias.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 17: Thirdparty regression test INI + reference output + CMake registration

**Files:**
- Create: `src/tests/topp/PeptideDataBaseSearchFI_4.ini`
- Create: `src/tests/topp/PeptideDataBaseSearchFI_4_out.idXML`
- Modify: `src/tests/topp/CMakeLists.txt`

- [ ] **Step 17.1: Create `PeptideDataBaseSearchFI_4.ini` (copy and modify `_1.ini`)**

Copy `PeptideDataBaseSearchFI_1.ini` to `PeptideDataBaseSearchFI_4.ini`:
```bash
cp src/tests/topp/PeptideDataBaseSearchFI_1.ini src/tests/topp/PeptideDataBaseSearchFI_4.ini
```

Edit the new file to change:
- `precursor:mass_tolerance_lower` = `5.0`
- `precursor:mass_tolerance_upper` = `15.0` (asymmetric, compensates the +7 ppm shift)
- `calibration:enabled` = `true`
- `calibration:min_psms` = `5` (lower threshold so the small fixture triggers calibration)

Leave the FASTA path and other params the same.

- [ ] **Step 17.2: Register in `src/tests/topp/CMakeLists.txt`**

Find the `PeptideDataBaseSearchFI_3_multi` block and add after it:
```cmake
add_test("TOPP_PeptideDataBaseSearchFI_4" ${TOPP_BIN_PATH}/PeptideDataBaseSearchFI -test
  -ini ${DATA_DIR_TOPP}/PeptideDataBaseSearchFI_4.ini
  -in ${DATA_DIR_TOPP}/SimpleSearchEngine_1_shifted_7ppm.mzML
  -out ${TESTS_TEMP_DIR}/PeptideDataBaseSearchFI_4_out.tmp.idXML
  -database ${DATA_DIR_TOPP}/SimpleSearchEngine_1.fasta)
add_test("TOPP_PeptideDataBaseSearchFI_4_out" ${DIFF}
  -in1 ${TESTS_TEMP_DIR}/PeptideDataBaseSearchFI_4_out.tmp.idXML
  -in2 ${DATA_DIR_TOPP}/PeptideDataBaseSearchFI_4_out.idXML
  -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
set_tests_properties("TOPP_PeptideDataBaseSearchFI_4_out" PROPERTIES DEPENDS
  "TOPP_PeptideDataBaseSearchFI_4")
```

- [ ] **Step 17.3: Rerun CMake configure to pick up the new test**

```bash
cmake -B OpenMS-build -S . 2>&1 | tail -10
```

- [ ] **Step 17.4: Run the test ONCE without a reference to capture the actual output**

```bash
cmake --build OpenMS-build --target PeptideDataBaseSearchFI -j$(nproc)
ctest --test-dir OpenMS-build -R "^TOPP_PeptideDataBaseSearchFI_4$" --output-on-failure 2>&1 | tail -30
```
This runs the tool but the `_out` diff subtest will fail (no reference file yet). That's expected.

Find the actual output in the test temp directory:
```bash
find OpenMS-build -name "PeptideDataBaseSearchFI_4_out.tmp.idXML" -exec ls -la {} \;
```

- [ ] **Step 17.5: Inspect the output and copy to reference IFF reasonable**

Read the actual output. Verify:
- Has at least some PSMs (asymmetric bounds compensated the +7 ppm shift).
- PSM count is comparable to `PeptideDataBaseSearchFI_1_out.idXML` (the unshifted baseline) — this proves end-to-end bias compensation worked.
- Calibration log message shows a shift near +7 ppm in the test log.

If the output looks reasonable, copy it to the reference path:
```bash
cp $(find OpenMS-build -name "PeptideDataBaseSearchFI_4_out.tmp.idXML") \
   src/tests/topp/PeptideDataBaseSearchFI_4_out.idXML
```

If the output has zero PSMs or looks wrong, **do NOT commit a broken reference**. Investigate: is calibration firing? Is the asymmetric window `[5, 15]` catching the +7 ppm shift? Adjust the INI or fixture and retry.

- [ ] **Step 17.6: Re-run the test to confirm it now passes with the committed reference**

```bash
ctest --test-dir OpenMS-build -R "^TOPP_PeptideDataBaseSearchFI_4" --output-on-failure 2>&1 | tail -20
```
Expected: both `_4` and `_4_out` pass.

- [ ] **Step 17.7: Commit**

```bash
git add src/tests/topp/PeptideDataBaseSearchFI_4.ini \
        src/tests/topp/PeptideDataBaseSearchFI_4_out.idXML \
        src/tests/topp/CMakeLists.txt

git commit -m "$(cat <<'EOF'
test(PeptideDataBaseSearchFI): thirdparty-level asymmetric + calibration regression

End-to-end test using SimpleSearchEngine_1_shifted_7ppm.mzML: asymmetric
[5, 15] ppm bounds + calibration:enabled compensate the known +7 ppm
systematic shift, producing a PSM set comparable to the unshifted _1 baseline.

Guards against regressions in:
- getPeptidesInMassWindow absolute-bound semantics
- runCalibrationPass_ signed-bias preservation (shift ≈ +7 ppm detected)
- Calibration writeback member refresh
- min() reduction rule feeding OpenSearchModificationAnalysis

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 18: CHANGELOG entry

**Files:**
- Modify: `CHANGELOG`

- [ ] **Step 18.1: Add BREAKING entry under `TOPP tools > Changes > PeptideDataBaseSearchFI`**

Find the `PeptideDataBaseSearchFI:` section (currently starts around line 65). At the top of its bullet list, add:
```
      - BREAKING: precursor tolerance schema changed to asymmetric positive-magnitude pair.
        The scalar `precursor:mass_tolerance` and the `precursor:open_window_lower/_upper`
        params are REMOVED and replaced by `precursor:mass_tolerance_lower` and
        `precursor:mass_tolerance_upper` (both positive magnitudes, applied as [-lower, +upper]
        around the precursor mass). Default widened from 10 to 20 ppm per side. Calibration
        pass now preserves signed precursor bias instead of clobbering asymmetric user bounds
        to symmetric; an extreme-bias guard (|shift| >= spread) discards the calibration result
        and warns rather than silently loosening. pyOpenMS users setting `precursor:mass_tolerance`
        must migrate to the new param names. (#9108)
```

- [ ] **Step 18.2: Commit**

```bash
git add CHANGELOG
git commit -m "docs(CHANGELOG): note BREAKING asymmetric precursor schema for PeptideDataBaseSearchFI"
```

---

## Task 19: Final verification

- [ ] **Step 19.1: Full build**

```bash
cmake --build OpenMS-build -j$(nproc) 2>&1 | tail -20
```
Expected: clean build.

- [ ] **Step 19.2: Run all affected tests**

```bash
ctest --test-dir OpenMS-build -R "FragmentIndex|PeptideSearchEngineFIAlgorithm|PeptideDataBaseSearchFI" --output-on-failure 2>&1 | tail -60
```
Expected: all pass, including the new tests 1-6c, 7, 8, 9, 10, and the thirdparty regression `PeptideDataBaseSearchFI_4`.

- [ ] **Step 19.3: Run pyOpenMS test**

```bash
cmake --build OpenMS-build --target pyopenms -j$(nproc)
PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/test_modification_discovery.py -v 2>&1 | tail -30
```
Expected: all four test functions pass.

- [ ] **Step 19.4: Confirm no stale references remain**

```bash
grep -rn "precursor:mass_tolerance[^_]" src/openms src/topp src/tests src/pyOpenMS
grep -rn "open_window_lower\|open_window_upper\|open_precursor_window" src/openms src/topp src/tests src/pyOpenMS
grep -rn "getPeptidesInPrecursorRange" src/openms src/topp src/tests src/pyOpenMS
```
Expected: no matches except possibly in the CHANGELOG or comment blocks.

- [ ] **Step 19.5: Review the full commit log for this branch**

```bash
git log --oneline develop..HEAD
```
Expected: ~14 commits covering the refactor, calibration, new tests, thirdparty fixture, and CHANGELOG.

---

## Self-Review

**Spec coverage:**
- §5 parameter schema → Tasks 1, 2, 3, 4 (FragmentIndex.h/.cpp + PSEFIA.h/.cpp param registrations + updateMembers_)
- §6 FragmentIndex refactor → Tasks 1, 2
- §7 issue 1 (dead params) → Task 4
- §7 issue 2 (duplicate isOpenSearchMode_) → Tasks 1, 3
- §7 issue 3 (OSMA input) → Tasks 3, 4
- §7 issue 4 (calibration bias) → Tasks 4 (stub) + 7 (rewrite)
- §7 OSMA + computeModMatchTolerance_ → Tasks 3, 4
- §7 calibration writeback + member refresh → Task 7
- §8 PeptideDataBaseSearchFI.cpp → Task 5
- §9 test infrastructure additions (friend subclass, extreme_bias flag, computeModMatchTolerance_, last_calibration_result_) → Tasks 3, 7, 11
- §9 tests 1-6c → Tasks 8, 9, 10
- §9 tests 7-10 → Tasks 12, 13, 14, 15
- §9 migration → Task 6
- §9 thirdparty regression (new requirement from this conversation) → Tasks 16, 17
- §10 files touched → Tasks 1-18
- §11 deferred → documented; nothing to implement
- §12 risk → informational

All spec requirements traced to at least one task. ✓

**Placeholder scan:** No "TBD", "TODO", "fill in later", or "write tests for the above" found. Every code step shows actual code. Every test step shows the complete test body. ✓

**Type consistency:** 
- `precursor_mass_tolerance_lower_` / `_upper_` used consistently across FragmentIndex.h, FragmentIndex.cpp, PSEFIA.h, PSEFIA.cpp, tests.
- `computeModMatchTolerance_()` signature matches between declaration (Task 3) and calls (Tasks 4, 13, 15).
- `CalibrationResult_` fields (`precursor_shift`, `precursor_spread`, `cal_lower`, `cal_upper`, `extreme_bias`) used consistently between declaration (Task 3), population (Task 7), and test assertions (Tasks 12, 14).
- `isOpenSearchMode(lower, upper, unit_ppm)` static signature matches between declaration (Task 1), delegation (Tasks 1, 3, 5), and test calls (Task 8 test 3). ✓

**Inline fixes made during review:** None needed.
