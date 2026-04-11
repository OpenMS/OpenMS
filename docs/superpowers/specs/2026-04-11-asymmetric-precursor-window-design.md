# Asymmetric Precursor Window for FragmentIndex / PeptideSearchEngineFI

**Date:** 2026-04-11
**Status:** Draft
**Tracks:** OpenMS/OpenMS#9108 item 5 (asymmetric precursor tolerance window)
**Related:** OpenMS/OpenMS#9107

## 1. Motivation

`PeptideSearchEngineFI` currently exposes a single symmetric `precursor:mass_tolerance` scalar (default 10 ppm) for closed search, plus a separate asymmetric `precursor:open_window_lower` / `precursor:open_window_upper` pair (default −100 / +200 Da) that only activates in open-search mode. The closed-search path hardcodes a `{0, 0}` window and uses `±mass_tolerance` symmetrically (`FragmentIndex.cpp:985-986`, `:812-813`).

Two independent forces push toward asymmetric closed-search tolerance:

1. **Comet / Sage / MSFragger parity.** All three reference engines expose asymmetric bounds as the primary tolerance knob.
2. **Wide-window acquisition (WWA) preparation.** WWA deposits fragments from a wide asymmetric isolation window into each MS2 spectrum; the search engine needs an asymmetric mass window that ultimately becomes per-spectrum (read from `Precursor::getIsolationWindow{Lower,Upper}Offset()`).

This spec unifies the two existing asymmetric-window mechanisms (closed-search tolerance + open-search window) into a single signed pair, refactors the underlying `FragmentIndex` API to an honest absolute-bounds interface, and declares a `wide_window` placeholder flag so the future per-spectrum path drops in without another refactor.

Because `PeptideSearchEngineFI` / `PeptideDataBaseSearchFI` have not yet shipped to users, we take the freedom to make breaking parameter changes rather than layer backward-compat shims.

## 2. Goals

- **G1.** Users can specify an asymmetric precursor window for closed search (e.g., `[−5, +15]` ppm to compensate a known calibration bias).
- **G2.** Closed and open search use one unified asymmetric window mechanism. No duplicate knobs.
- **G3.** The `FragmentIndex` candidate-lookup API has a single, honest, absolute-bounds signature. The current additive `prec_tol + window` semantics — which hide the tolerance inside the function and only accept an *additional* offset from the caller — are removed.
- **G4.** `wide_window` flag declared (default `false`, implementation deferred) so the per-spectrum WWA path can wire into the same refactored interface without further churn.
- **G5.** No regression in modification discovery for today's open-search flow.

## 3. Non-goals

- **NG1.** Implementing wide-window acquisition / per-spectrum isolation-window candidate selection. The flag is declared; the implementation is a separate design.
- **NG2.** Chimeric-spectrum support (issue #9108 item 16). Separate design.
- **NG3.** A dedicated, tight tolerance for `OpenSearchModificationAnalysis` mod-delta → UniMod matching (see §7). Preserved as-is; future refinement.
- **NG4.** Backward compatibility with the existing `precursor:mass_tolerance` scalar or `precursor:open_window_*` params. Clean removal. Test INIs are migrated.

## 4. Reference: how other engines handle this

| Engine | Schema | Signed? | Open/wide handling |
|---|---|---|---|
| **Sage** | `precursor_tol: { ppm: [lo, hi] }` or `{ da: [lo, hi] }` | yes | `wide_window: true` ignores `precursor_tol` and reads per-spectrum isolation windows |
| **MSFragger** | `precursor_mass_lower` / `precursor_mass_upper` (Da) | yes (lower ∈ [−150, 0]) | Same pair scales up; MSFragger-DDA+ adds per-spectrum WWA |
| **Comet 2024+** | Legacy `peptide_mass_tolerance` **plus** `peptide_mass_tolerance_lower` / `_upper` | positive magnitudes | Separate knobs for open-mode behavior |

Sage and MSFragger share the clean model: **one signed pair, one mechanism, wide-window is an explicit flag that switches the bounds source**. We follow this.

## 5. Parameter schema

### Added

| Param | Type | Default | Notes |
|---|---|---|---|
| `precursor:mass_tolerance_lower` | `double` | `−10.0` | **Signed.** The lower offset from the precursor mass. Must be ≤ 0. Unit controlled by `precursor:mass_tolerance_unit`. |
| `precursor:mass_tolerance_upper` | `double` | `+10.0` | **Signed.** The upper offset from the precursor mass. Must be ≥ 0. Same unit. |
| `precursor:wide_window` | `string` (enum `"true"`/`"false"`) | `"false"` | **Declared but not implemented in this PR.** When `true`, the search ignores `mass_tolerance_lower`/`_upper` and reads per-spectrum isolation-window offsets. Setting `true` logs a warning and throws `NotImplementedError`. |

### Removed

| Param | Replacement |
|---|---|
| `precursor:mass_tolerance` | `mass_tolerance_lower` + `mass_tolerance_upper` |
| `precursor:open_window_lower` | `mass_tolerance_lower` |
| `precursor:open_window_upper` | `mass_tolerance_upper` |

### Retained

| Param | Notes |
|---|---|
| `precursor:mass_tolerance_unit` | Still `"ppm"` | `"Da"`. Applies to both signed bounds. |
| `precursor:isotope_error_min` / `_max` | Unchanged. Still iterated in **closed** mode only (§6 — `searchDifferentPrecursorRanges`). |
| `precursor:isotopes` (PeptideSearchEngineFIAlgorithm) | Unchanged. Pre-search peak-misassignment correction, independent of the tolerance window. |

### Validation

- `mass_tolerance_lower ≤ 0 ≤ mass_tolerance_upper`. A nonzero-width window (`upper − lower > 0`) is required.
- Unit `"Da"` accepts any magnitude; unit `"ppm"` is unbounded but auto-trips open-mode detection (§6) above threshold.

### Open-search auto-detection

`isOpenSearchMode_()` returns `true` iff
```
max(|lower|, upper) > threshold
```
where `threshold = 1000.0` (ppm unit) or `1.0` (Da unit). Same thresholds as today, applied to the effective asymmetric bounds.

## 6. Interface refactor: `FragmentIndex`

### Current signature

```cpp
// FragmentIndex.h:197
std::pair<size_t, size_t> getPeptidesInPrecursorRange(float precursor_mass,
                                                     const std::pair<float, float>& window);
```
Currently computes
```
left  = precursor_mass − prec_tol + window.first
right = precursor_mass + prec_tol + window.second
```
(`FragmentIndex.cpp:810-814`). The `window` is an *additive offset* on top of the implicit symmetric `±prec_tol`. Callers must know this. The name does not reveal it.

### New signature

```cpp
// FragmentIndex.h
/// Return the [begin, end) peptide index range whose precursor mass falls in
/// [mass + window.first, mass + window.second]. Bounds are absolute signed
/// offsets — no hidden tolerance is added. window.first must be ≤ 0, window.second must be ≥ 0.
std::pair<size_t, size_t> getPeptidesInMassWindow(float precursor_mass,
                                                  const std::pair<float, float>& window) const;
```

Implementation:
```cpp
auto left_it  = std::lower_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                 precursor_mass + window.first,
                                 [](const Peptide& a, float b) { return a.precursor_mz_ < b; });
auto right_it = std::upper_bound(fi_peptides_.begin(), fi_peptides_.end(),
                                 precursor_mass + window.second,
                                 [](float b, const Peptide& a) { return b < a.precursor_mz_; });
```
No `prec_tol` computation inside the function. Callers own the unit conversion.

### New members (`FragmentIndex.h`)

Replace:
```cpp
float precursor_mz_tolerance_;
bool  precursor_mz_tolerance_unit_ppm_{true};
double open_precursor_window_lower_;
double open_precursor_window_upper_;
```
with:
```cpp
double precursor_mass_tolerance_lower_{-10.0};   // signed, ≤ 0
double precursor_mass_tolerance_upper_{+10.0};   // signed, ≥ 0
bool   precursor_mass_tolerance_unit_ppm_{true};
bool   wide_window_{false};                      // placeholder
```

### New helper

```cpp
/// Compute the absolute mass window [lo, hi] around a given precursor_mass,
/// converting ppm → Da at that mass if the unit is ppm. Returns signed offsets.
std::pair<float, float> computeMassWindow_(float precursor_mass) const;
```
Body:
```cpp
if (precursor_mass_tolerance_unit_ppm_) {
  float lo = Math::ppmToMass((float)precursor_mass_tolerance_lower_, precursor_mass);
  float hi = Math::ppmToMass((float)precursor_mass_tolerance_upper_, precursor_mass);
  return {lo, hi};   // lo is negative because lower_ is negative
}
return {(float)precursor_mass_tolerance_lower_, (float)precursor_mass_tolerance_upper_};
```

### `searchDifferentPrecursorRanges` rewrite

```cpp
void FragmentIndex::searchDifferentPrecursorRanges(const MSSpectrum& spectrum,
                                                   float precursor_mass,
                                                   SpectrumMatchesTopN& sms,
                                                   uint16_t charge)
{
  // Open mode absorbs isotope shifts into the wide window (same behavior as today).
  const bool open_mode = isOpenSearchMode_();
  const int16_t iso_lo = open_mode ? 0 : min_isotope_error_;
  const int16_t iso_hi = open_mode ? 0 : max_isotope_error_;

  for (int16_t isotope_error = iso_lo; isotope_error <= iso_hi; ++isotope_error)
  {
    const float shifted_mass = precursor_mass
      + (float)isotope_error * (float)Constants::C13C12_MASSDIFF_U;

    const auto window = computeMassWindow_(shifted_mass);

    SpectrumMatchesTopN candidates;
    auto range = getPeptidesInMassWindow(shifted_mass, window);
    candidates.hits_.resize(range.second - range.first + 1);
    queryPeaks(candidates, spectrum, range, isotope_error, charge);
    sms += candidates;
  }
}
```

No additive mixing of tolerance and window. No `{0, 0}` hardcode. One code path for closed and open modes. Isotope-error iteration-to-[0,0] behavior in open mode is preserved.

### `updateMembers_()` changes

Replace the block at `FragmentIndex.cpp:1194-1219`:

```cpp
precursor_mass_tolerance_lower_     = param_.getValue("precursor:mass_tolerance_lower");
precursor_mass_tolerance_upper_     = param_.getValue("precursor:mass_tolerance_upper");
precursor_mass_tolerance_unit_ppm_  = param_.getValue("precursor:mass_tolerance_unit").toString() == "ppm";
wide_window_                        = param_.getValue("precursor:wide_window").toString() == "true";

// Validation
if (precursor_mass_tolerance_lower_ > 0.0 || precursor_mass_tolerance_upper_ < 0.0)
{
  throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    "precursor:mass_tolerance_lower must be ≤ 0 and mass_tolerance_upper ≥ 0");
}
if (wide_window_)
{
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
}

if (isOpenSearchMode_())
{
  OPENMS_LOG_INFO << "[FragmentIndex] Open-search mode (window ["
                  << precursor_mass_tolerance_lower_ << ", "
                  << precursor_mass_tolerance_upper_ << "] "
                  << (precursor_mass_tolerance_unit_ppm_ ? "ppm" : "Da") << ")." << std::endl;
}
```

### `isOpenSearchMode_()` update

```cpp
// FragmentIndex.h (inline)
bool isOpenSearchMode_() const
{
  const double threshold = precursor_mass_tolerance_unit_ppm_ ? 1000.0 : 1.0;
  return std::max(std::abs(precursor_mass_tolerance_lower_),
                  precursor_mass_tolerance_upper_) > threshold;
}
```

## 7. Call-site changes: `PeptideSearchEngineFIAlgorithm`

Three issues found in the existing code that this refactor cleans up:

1. **Dead parameter registration** (`PeptideSearchEngineFIAlgorithm.cpp:181-182`): registers `precursor:open_window_lower`/`_upper` but never reads them in `updateMembers_()`. Pure pass-through to `FragmentIndex`. **Delete these two lines.**
2. **Duplicate `isOpenSearchMode_()`** (`PeptideSearchEngineFIAlgorithm.h:478-483`): a second copy of the auto-detection logic lives as an inline method in the header, alongside `FragmentIndex`'s. It currently reads `precursor_mass_tolerance_` (scalar) and `precursor_mass_tolerance_unit_`. Update to the new asymmetric-bounds rule.
3. **`OpenSearchModificationAnalysis` tolerance input** (`PeptideSearchEngineFIAlgorithm.cpp:889-895`): currently passes `precursor_mass_tolerance_` (the scalar) and the unit. It uses this tolerance to match observed Δmass to UniMod entries.

### Parameter registration changes (PeptideSearchEngineFIAlgorithm.cpp:65-72, 181-182)

Replace
```cpp
defaults_.setValue("precursor:mass_tolerance", 10.0, "+/- tolerance for precursor mass.");
```
with
```cpp
defaults_.setValue("precursor:mass_tolerance_lower", -10.0,
                   "Lower (signed, ≤ 0) offset from precursor mass for candidate selection.");
defaults_.setValue("precursor:mass_tolerance_upper", +10.0,
                   "Upper (signed, ≥ 0) offset from precursor mass for candidate selection.");
defaults_.setValue("precursor:wide_window", "false",
                   "Declared for forward compatibility — per-spectrum isolation window search. "
                   "NOT YET IMPLEMENTED.");
defaults_.setValidStrings("precursor:wide_window", {"true", "false"});
```
Delete lines 181-182 (`open_window_lower` / `_upper` registrations — dead).

### `updateMembers_()` changes (cpp:223)

Replace
```cpp
precursor_mass_tolerance_ = param_.getValue("precursor:mass_tolerance");
```
with
```cpp
precursor_mass_tolerance_lower_ = param_.getValue("precursor:mass_tolerance_lower");
precursor_mass_tolerance_upper_ = param_.getValue("precursor:mass_tolerance_upper");
```

### `isOpenSearchMode_()` duplicate update (`PeptideSearchEngineFIAlgorithm.h:478-483`)

Replace the inline body:
```cpp
bool isOpenSearchMode_() const
{
  const double threshold = precursor_mass_tolerance_unit_ == "ppm" ? 1000.0 : 1.0;
  return std::max(std::abs(precursor_mass_tolerance_lower_),
                  precursor_mass_tolerance_upper_) > threshold;
}
```
This duplicates `FragmentIndex::isOpenSearchMode_()` by necessity — both classes have their own param copies. Keeping them mechanically identical is the minimum bar; consolidating into a shared helper is listed as a follow-up (§11).

### Calibration pass changes (cpp:760-782)

The calibration pass currently writes a scalar back:
```cpp
fi_params.setValue("precursor:mass_tolerance", effective_precursor_tol);
```
New behavior: write back symmetric bounds derived from the calibrated scalar:
```cpp
fi_params.setValue("precursor:mass_tolerance_lower", -effective_precursor_tol);
fi_params.setValue("precursor:mass_tolerance_upper", +effective_precursor_tol);
```
The calibration result is (and stays) a single scalar — calibration estimates symmetric instrument precision, not asymmetric bias. If a user configured asymmetric bounds and calibration succeeds, the calibration result **overrides** to symmetric. This is the correct behavior: calibration finds the instrument's precision, which the user's asymmetric-bias setting cannot improve on. A log line should explain the override.

### `OpenSearchModificationAnalysis` input (cpp:889-895)

Change from
```cpp
mod_analyzer.analyzeModifications(
  peptide_ids,
  precursor_mass_tolerance_,
  precursor_mass_tolerance_unit_ == "ppm",
  false,
  "");
```
to
```cpp
// OpenSearchModificationAnalysis matches observed Δmass to UniMod using a scalar
// tolerance. Under asymmetric bounds there is no single "the tolerance" value; we
// feed the wider side (conservative — same as today for the symmetric case, no
// regression for current open-search flows). A dedicated mod-matching tolerance
// (tighter than the search window) is deferred future work.
const double mod_match_tol = std::max(std::abs(precursor_mass_tolerance_lower_),
                                      precursor_mass_tolerance_upper_);
mod_analyzer.analyzeModifications(
  peptide_ids,
  mod_match_tol,
  precursor_mass_tolerance_unit_ == "ppm",
  false,
  "");
```

No change to `OpenSearchModificationAnalysis` itself. The semantic-regression risk is bounded: for symmetric configurations the value is identical; for today's open-search flow with `mass_tolerance > 1 Da`, the user's existing scalar becomes `max(|lower|, upper)` — which is the same number, by construction.

## 8. Call-site changes: `PeptideDataBaseSearchFI.cpp`

The TOPP tool wrapper reads `precursor:mass_tolerance` and `mass_tolerance_unit` directly at `PeptideDataBaseSearchFI.cpp:168-169` and recomputes open-search mode locally at `:170-172`. Update both to read the new pair and apply the new detection rule. No other touches — the tool copies the full `Search:` param subtree to the algorithm and lets the algorithm handle the rest.

## 9. Tests

### Migration (mechanical)

- `src/tests/topp/PeptideDataBaseSearchFI_1.ini`: `mass_tolerance = 5` → `mass_tolerance_lower = -5`, `mass_tolerance_upper = +5`.
- `src/tests/topp/PeptideDataBaseSearchFI_2.ini`: `mass_tolerance = 10.0` → `mass_tolerance_lower = -10.0`, `mass_tolerance_upper = +10.0`.
- `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`: every site that sets `precursor:mass_tolerance` (grep at implementation time) → replace with the signed pair. Output references should be **unchanged** — if they differ, that's a real regression to investigate, not a reference to update.

### New tests

1. **Symmetric default equivalence**: build the index at the default `[−10, +10]` ppm, search a fixture spectrum, confirm the PSM set is identical to the pre-refactor baseline. (Guards the refactor.)
2. **Asymmetric closed**: build the index with `[−5, +15]` ppm on a synthetic peptide mix where the precursor masses are shifted by `+8 ppm`. Confirm identification succeeds with the asymmetric bounds and fails with `[−5, +5]` ppm (calibration-offset simulation).
3. **Open-mode auto-detection boundary**: confirm `[−500, +1500]` ppm auto-flips to open mode; `[−999, +999]` ppm stays closed.
4. **Open-mode isotope-error force-to-zero**: confirm the iteration still collapses to `[0, 0]` in open mode (same behavior as today).
5. **Validation**: `mass_tolerance_lower = +5` throws `InvalidParameter`.
6. **`wide_window = "true"` throws `NotImplemented`**.

### Regression watch

Run the full `FragmentIndex` and `PeptideSearchEngineFIAlgorithm` class tests + the two TOPP `PeptideDataBaseSearchFI_*` integration tests after the refactor. Any diff in the idXML reference outputs for the symmetric tests must be investigated as a real regression, not bulk-accepted.

## 10. Files touched

| File | Change |
|---|---|
| `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h` | rename `getPeptidesInPrecursorRange` → `getPeptidesInMassWindow`; add `computeMassWindow_`; replace members; update `isOpenSearchMode_()` |
| `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` | implementations; rewrite `searchDifferentPrecursorRanges`; `updateMembers_()`; param registrations; delete dead `modification window` comment |
| `src/openms/include/OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h` | replace scalar tolerance members with `_lower_` / `_upper_` |
| `src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp` | param registrations (add pair, delete dead `open_window_*`); `updateMembers_()`; local `isOpenSearchMode_()`; calibration pass writeback; `OpenSearchModificationAnalysis` tolerance input |
| `src/topp/PeptideDataBaseSearchFI.cpp` | param reads + local open-mode detection |
| `src/tests/topp/PeptideDataBaseSearchFI_1.ini` / `_2.ini` | migrate |
| `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp` | migrate + add new tests |
| `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` | add tests (2)–(6) above |
| `CHANGELOG` | breaking change entry |

## 11. Deferred / follow-ups

- **Wide-window implementation**: per-spectrum `Precursor::getIsolationWindow{Lower,Upper}Offset()` → `getPeptidesInMassWindow` call path. Separate design, but the shape of `getPeptidesInMassWindow` is already correct.
- **Dedicated mod-matching tolerance**: `modification:matching_tolerance` param for `OpenSearchModificationAnalysis`, decoupled from the (wide) search window. Improves UniMod matching precision in open mode. Separate design.
- **Consolidating the duplicate `isOpenSearchMode_()`** between `FragmentIndex` and `PeptideSearchEngineFIAlgorithm` into a single helper. Minor cleanup.

## 12. Risk

- **Low** — changes are confined to one search engine with no released users. The `FragmentIndex` API rename affects only `PeptideSearchEngineFIAlgorithm` and the class test; grepped, no other callers.
- **Test reference drift is the main way this can silently regress.** The mitigation is: symmetric-default test must show no idXML output diff; any diff is treated as a real regression and investigated root-cause.
