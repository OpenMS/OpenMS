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

This spec unifies the two existing asymmetric-window mechanisms (closed-search tolerance + open-search window) into a single positive-magnitude pair and refactors the underlying `FragmentIndex` API to an honest absolute-bounds interface. Wide-window acquisition (WWA) is listed as follow-up work in §11; **no `wide_window` flag is declared in this PR** — it would be premature commitment to a param name before the WWA design is nailed down, and adds param-schema churn for zero current user value.

Because `PeptideSearchEngineFI` / `PeptideDataBaseSearchFI` have not yet shipped to users, we take the freedom to make breaking parameter changes rather than layer backward-compat shims.

## 2. Goals

- **G1.** Users can specify an asymmetric precursor window for closed search (e.g., `mass_tolerance_lower=5, mass_tolerance_upper=15` ppm → effective window `[−5, +15]` ppm to compensate a known calibration bias). Calibration pass (§7) preserves this bias instead of clobbering it to symmetric.
- **G2.** Closed and open search use one unified asymmetric window mechanism. No duplicate knobs.
- **G3.** The `FragmentIndex` candidate-lookup API has a single, honest, absolute-bounds signature. The current additive `prec_tol + window` semantics — which hide the tolerance inside the function and only accept an *additional* offset from the caller — are removed.
- **G4.** No regression in modification discovery for today's open-search flow.
- **G5.** Calibration pass preserves asymmetric user bounds — the flagship use case (a user compensating a known instrument bias via `[15, 5]` ppm) must not have its window clobbered back to symmetric by the calibration step.

## 3. Non-goals

- **NG1.** Wide-window acquisition (WWA) in any form. No `wide_window` param, no per-spectrum isolation-window candidate selection, no placeholder flag. The WWA work is listed in §11 as a dedicated follow-up design; committing to a flag name today before the implementation shape is settled would be premature.
- **NG2.** Chimeric-spectrum support (issue #9108 item 16). Separate design.
- **NG3.** A dedicated, tight tolerance for `OpenSearchModificationAnalysis` mod-delta → UniMod matching. Preserved as-is; `OpenSearchModificationAnalysis` currently clamps every input at `MAX_MOD_MAPPING_TOL_ = 0.02 Da` (and discards ppm inputs entirely), so the reduction rule is behavior-neutral today (§7). Future work listed in §11.
- **NG4.** Backward compatibility with the existing `precursor:mass_tolerance` scalar or `precursor:open_window_*` params. Clean removal. Test INIs are migrated.

## 4. Reference: how other engines handle this

| Engine | Schema | Sign convention | Open/wide handling |
|---|---|---|---|
| **Sage** | `precursor_tol: { ppm: [lo, hi] }` or `{ da: [lo, hi] }` | signed (`lo < 0 < hi`) | `wide_window: true` ignores `precursor_tol` and reads per-spectrum isolation windows |
| **MSFragger (raw config)** | `precursor_mass_lower` / `precursor_mass_upper` (Da) | signed (`lower ∈ [−150, 0]`) | Same pair scales up; MSFragger-DDA+ adds per-spectrum WWA |
| **OpenMS `MSFraggerAdapter`** | `precursor_mass_tolerance_lower` / `_upper` | **positive magnitudes**, sign flipped at write-time (`-arg_precursor_mass_tolerance_lower`) | Same pair scales up |
| **Comet 2024+** | Legacy `peptide_mass_tolerance` **plus** `peptide_mass_tolerance_lower` / `_upper` | positive magnitudes | Separate knobs for open-mode behavior |

**Convention chosen: positive magnitudes**, matching `MSFraggerAdapter.cpp:97-98,252-253,987-988` and Comet. The internal FragmentIndex helper applies the signs (`[-lower, +upper]`). Rationale:

- Keeps the user-facing config form consistent across `MSFraggerAdapter`, `CometAdapter`, and `PeptideDataBaseSearchFI`. Two OpenMS adapters disagreeing on the sign convention for identically named params would be a real footgun.
- Users rarely need truly asymmetric windows in practice; magnitudes read naturally as "how many ppm/Da of tolerance on each side".
- Sign flips happen exactly once, inside `computeMassWindow_()`, where the convention is centralized.

## 5. Parameter schema

### Added

| Param | Type | Default | Notes |
|---|---|---|---|
| `precursor:mass_tolerance_lower` | `double` | `20.0` | **Positive magnitude.** Effective lower bound is `−lower`. Must be finite and ≥ 0. Unit controlled by `precursor:mass_tolerance_unit`. Default widened from the original `10.0` to reduce first-run surprises on uncalibrated / non-Orbitrap instruments (Bruker timsTOF typically needs 15–20 ppm). The internal calibration pass (§7) tightens it on real data. |
| `precursor:mass_tolerance_upper` | `double` | `20.0` | **Positive magnitude.** Effective upper bound is `+upper`. Must be finite and ≥ 0. Same unit. |

### Removed

| Param | Replacement |
|---|---|
| `precursor:mass_tolerance` | `mass_tolerance_lower` + `mass_tolerance_upper` |
| `precursor:open_window_lower` | `mass_tolerance_lower` |
| `precursor:open_window_upper` | `mass_tolerance_upper` |

### Retained

| Param | Notes |
|---|---|
| `precursor:mass_tolerance_unit` | Unchanged enum (`"ppm"` or `"Da"`). Applies to both magnitudes. |
| `precursor:isotope_error_min` / `_max` | Unchanged. Still iterated in **closed** mode only (§6 — `searchDifferentPrecursorRanges`). |
| `precursor:isotopes` (PeptideSearchEngineFIAlgorithm) | Unchanged. Pre-search peak-misassignment correction, independent of the tolerance window. |

### Validation

- Both `mass_tolerance_lower` and `mass_tolerance_upper` must be **finite** and `>= 0` (positive magnitudes). `setMinFloat(0.0)` at registration catches negatives via `Param::checkDefaults_`, but does **not** reject `+inf` or `NaN` (see `ParamEntry::isValid` at `Param.cpp:143-166` — `NaN < 0.0` is `false`, so NaN passes the min check). `FragmentIndex::updateMembers_()` therefore adds an explicit `std::isfinite` guard — NaN would otherwise reach `lower_bound`'s comparator and break strict-weak-ordering.
- Nonzero total width (`lower + upper > 0`) is required — checked in `updateMembers_()` after default merging, since `setMinFloat(0.0)` accepts zero.
- Unit `"Da"` accepts any magnitude; unit `"ppm"` is unbounded but auto-trips open-mode detection below.
- Zero on one side is legal: `[0, 5]` ppm means "the target mass plus 0 to 5 ppm" — a one-sided half-line window, useful when the user wants to encode "only search *above* the observed precursor".

### Open-search auto-detection

`isOpenSearchMode_()` returns `true` iff
```
max(lower, upper) > threshold
```
where `threshold = 1000.0` (ppm unit) or `1.0` (Da unit). Same thresholds as today, applied to the positive magnitudes. The check is strict `>`, so exactly `1000.0` ppm stays closed-mode (preserves today's behavior).

**Warn on silent mode flip.** When a user ports a wide config (e.g., `[2000, 2000]` ppm meant as loose closed search) that trips the threshold, `updateMembers_()` emits `OPENMS_LOG_WARN` explicitly naming the mode switch and the downstream consequence (isotope-error iteration collapses to `[0, 0]`). A plain INFO line is easy to miss.

### Isotope-error interaction note

Under asymmetric closed bounds (e.g., `lower=15, upper=5` ppm for a negative calibration bias), the default `precursor:isotope_error_min = 0` silently misses candidates where the monoisotopic peak was mis-picked **upward** one `13C` — a case a symmetric 15 ppm tolerance with isotope iteration would catch. **Users configuring strongly asymmetric closed bounds should also review `isotope_error_min`.** The param default is not changed in this PR (out of scope), but the param description is updated to call out the interaction.

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
/// Return the [begin_idx, end_idx) range such that
/// fi_peptides_[i].precursor_mz_ ∈ [mass + window.first, mass + window.second]
/// for all i ∈ [begin_idx, end_idx).
/// @p window holds **signed absolute offsets** (the first element is ≤ 0, the second ≥ 0).
/// No hidden tolerance is added — the caller supplies the full signed window.
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
Unit conversion lives in `computeMassWindow_`, not in the bounds accessor.

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
double precursor_mass_tolerance_lower_{20.0};    // positive magnitude, applied as -lower
double precursor_mass_tolerance_upper_{20.0};    // positive magnitude, applied as +upper
bool   precursor_mass_tolerance_unit_ppm_{true};
```

### Includes

`FragmentIndex.h` currently lacks `<algorithm>`. The new inline static `isOpenSearchMode` uses `std::max`. Add:
```cpp
#include <algorithm>   // std::max
```
`PeptideSearchEngineFIAlgorithm.h` only needs to include `FragmentIndex.h` (which it already does transitively) to reach the static helper. The `.cpp` files that call `std::abs` for the calibration-pass math already include `<cmath>` via existing paths — verify during implementation.

### New helper

```cpp
/// Compute the signed mass window {lo, hi} around a given precursor_mass,
/// converting ppm → Da at that mass if the unit is ppm. lo is negative, hi is positive.
std::pair<float, float> computeMassWindow_(float precursor_mass) const;
```
Body:
```cpp
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
```

Notes for implementers:
- Explicit `Math::ppmToMass<float>(...)` — the template deduces a single `T` from both args, so mixing `float` and `double` fails to compile. Pin to `float` since `precursor_mass` and the stored peptide masses are `float`.
- `static_cast<float>` over C-style casts (OpenMS convention).
- Sign application is centralized here — this is the only place positive magnitudes become signed bounds.

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
      + static_cast<float>(isotope_error) * static_cast<float>(Constants::C13C12_MASSDIFF_U);

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

// Validation — `setMinFloat(0.0)` rejects negatives via checkDefaults_, but NaN and +inf
// slip through (NaN < 0 is false). NaN would break lower_bound's strict-weak-ordering.
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

### `isOpenSearchMode_()` update + deduplication

The existing duplicate implementations in `FragmentIndex` and `PeptideSearchEngineFIAlgorithm` — root cause of bug #2 in §7 — are consolidated into one inline static helper on `FragmentIndex`, then both classes call it.

```cpp
// FragmentIndex.h (inline static)
static bool isOpenSearchMode(double lower_magnitude,
                             double upper_magnitude,
                             bool unit_ppm) noexcept
{
  const double threshold = unit_ppm ? 1000.0 : 1.0;
  return std::max(lower_magnitude, upper_magnitude) > threshold;
}

// Convenience non-static wrapper used inside FragmentIndex itself:
bool isOpenSearchMode_() const noexcept
{
  return isOpenSearchMode(precursor_mass_tolerance_lower_,
                          precursor_mass_tolerance_upper_,
                          precursor_mass_tolerance_unit_ppm_);
}
```

`PeptideSearchEngineFIAlgorithm::isOpenSearchMode_()` (header inline) calls the static helper with its own members (§7). Single source of truth, no `std::abs` (magnitudes are already positive), follow-up §11 item removed.

### `const`-correctness chain

`getPeptidesInMassWindow`, `computeMassWindow_`, `isOpenSearchMode_()` are all marked `const`. Without this, the public `const` method can't call the private helpers. This is a minimal fix — the broader `queryPeaks` / `querySpectrum` / `searchDifferentPrecursorRanges` chain remains non-`const` (out of scope for this refactor).

## 7. Call-site changes: `PeptideSearchEngineFIAlgorithm`

Four issues found in the existing code that this refactor cleans up:

1. **Dead parameter registration** (`PeptideSearchEngineFIAlgorithm.cpp:181-182`): registers `precursor:open_window_lower`/`_upper` but never reads them in `updateMembers_()`. Pure pass-through to `FragmentIndex`. Delete.
2. **Duplicate `isOpenSearchMode_()`** (`PeptideSearchEngineFIAlgorithm.h:478-483`): a second copy of the auto-detection logic lives as an inline method in the header, alongside `FragmentIndex`'s. Consolidate via the new static helper on `FragmentIndex`.
3. **`OpenSearchModificationAnalysis` tolerance input** (`PeptideSearchEngineFIAlgorithm.cpp:889-895`): currently passes `precursor_mass_tolerance_` (scalar) and unit. The scalar is used to match observed Δmass to UniMod entries.
4. **Calibration pass discards signed bias** (`PeptideSearchEngineFIAlgorithm.cpp:1580`): `runCalibrationPass_` computes signed `prec_err` at `:1534` then throws the sign away via `std::abs` at `:1580`. Writeback at `:770-780` updates **only** the `fragment_index_` params, leaving the algo-level `precursor_mass_tolerance_*_` stale — a latent bug today (since `OpenSearchModificationAnalysis` reads the stale members), which would silently survive the refactor.

### Parameter registration changes (PeptideSearchEngineFIAlgorithm.cpp:65-72, 181-182)

Replace
```cpp
defaults_.setValue("precursor:mass_tolerance", 10.0, "+/- tolerance for precursor mass.");
```
with
```cpp
defaults_.setValue("precursor:mass_tolerance_lower", 20.0,
                   "Lower-side precursor-mass tolerance (positive magnitude; effective window "
                   "is [-lower, +upper] around the precursor). "
                   "When strongly asymmetric, also review precursor:isotope_error_min.");
defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
defaults_.setValue("precursor:mass_tolerance_upper", 20.0,
                   "Upper-side precursor-mass tolerance (positive magnitude).");
defaults_.setMinFloat("precursor:mass_tolerance_upper", 0.0);
```
Delete lines 181-182 (`open_window_lower` / `_upper` registrations — dead).

Section description at `PeptideSearchEngineFIAlgorithm.cpp:77` should mention the positive-magnitude convention explicitly so users don't enter negative numbers:
```cpp
defaults_.setSectionDescription("precursor",
  "Precursor (Parent Ion) Options. mass_tolerance_lower/_upper are positive magnitudes "
  "applied as [-lower, +upper] around the precursor mass.");
```

**Mirror the same registrations in `FragmentIndex::defaults_`** at `FragmentIndex.cpp:1083-1088`. The dual registration between `FragmentIndex` and `PeptideSearchEngineFIAlgorithm` is the existing OpenMS pattern — the algorithm copies its full `Search:` param subtree to `FragmentIndex::setParameters()` at runtime (`PeptideSearchEngineFIAlgorithm.cpp:695`), so both sides must declare the same defaults to satisfy `checkDefaults_()`. This refactor preserves that pattern and does not try to consolidate.

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

The finiteness / zero-width validation lives in `FragmentIndex::updateMembers_()` (§6) and fires when the algorithm's param subtree is pushed into `fragment_index_.setParameters()` at `cpp:695`. No duplicate validation needed here.

### `isOpenSearchMode_()` duplicate → shared helper (`PeptideSearchEngineFIAlgorithm.h:478-483`)

Replace the inline body with a call to the shared static helper on `FragmentIndex`:
```cpp
bool isOpenSearchMode_() const
{
  return FragmentIndex::isOpenSearchMode(precursor_mass_tolerance_lower_,
                                          precursor_mass_tolerance_upper_,
                                          precursor_mass_tolerance_unit_ == "ppm");
}
```
One source of truth. This was flagged in the spec review as "cheap to fix now" rather than deferred — it's the root cause of the two-classes-drifting risk and adds no header coupling beyond what's already there (`PeptideSearchEngineFIAlgorithm.h` already includes `FragmentIndex.h`).

### Calibration pass — preserve signed bias (cpp:760-782, 1454-1626)

The current `runCalibrationPass_` computes signed `prec_err` at cpp:1534 then throws the sign away via `std::abs` at cpp:1580 and returns a single-scalar spread. Writing `[−spread, +spread]` back discards any real calibration bias — defeating the whole point of asymmetric closed search. A user setting `[15, 5]` ppm to compensate a −5 ppm bias would have the window clobbered to `[3, 3]` ppm and lose every target.

**Fix summary.** Extend `CalibrationResult_` with a signed shift plus a spread-around-the-shift and an `extreme_bias` flag; compute them from the signed residuals; write asymmetric bounds back **unless** `|shift| > spread` (extreme bias, positive-magnitude schema can't faithfully represent a one-sided window without a lossy approximation that would widen the result — we skip instead and log a WARN).

Extended `CalibrationResult_` (`PeptideSearchEngineFIAlgorithm.h:444-450`):
```cpp
struct CalibrationResult_
{
  double precursor_shift{0};     ///< median of signed precursor errors (bias)
  double precursor_spread{0};    ///< median(|e - shift|) + 3*MAD(|e - shift|)
  double cal_lower{0};           ///< calibrated lower magnitude (valid iff !extreme_bias)
  double cal_upper{0};           ///< calibrated upper magnitude (valid iff !extreme_bias)
  double fragment_tolerance{0};
  double fragment_shift{0};
  bool   extreme_bias{false};    ///< |shift| > spread — writeback skipped (test hook)
  bool   success{false};
};
```
`precursor_tolerance` is replaced. Verified by the header-coupling review: nothing outside `PeptideSearchEngineFIAlgorithm.{h,cpp}` references `CalibrationResult_::precursor_tolerance` (no pyOpenMS bindings, no tests). Safe to remove.

Updated `runCalibrationPass_` (replacing cpp:1575-1601). Uses the existing `min_tolerance` local at `cpp:1582` — no new constant. The formula `res_median + 3*MAD(residuals)` mirrors the current `median(|e|) + 3*MAD(|e|)` structure applied to residuals around the signed shift rather than raw errors, so the symmetric-case spread is unchanged for a zero-centered distribution:

```cpp
// Signed median gives the calibration bias; spread is median(|e - shift|) + 3*MAD(|e - shift|),
// i.e., the same "typical residual + 3*scale" shape as the pre-refactor formula,
// just applied to residuals around the signed center instead of raw errors.
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
  // Only tighten. True calibrated signed window is [shift - spread, shift + spread],
  // which straddles zero because |shift| < spread. Map to positive magnitudes:
  //   -cal_lower = shift - spread  =>  cal_lower = spread - shift (> 0)
  //   +cal_upper = shift + spread  =>  cal_upper = spread + shift (> 0)
  const double cal_lower_raw = result.precursor_spread - prec_shift;
  const double cal_upper_raw = result.precursor_spread + prec_shift;
  result.cal_lower = std::min(cal_lower_raw, precursor_mass_tolerance_lower_);
  result.cal_upper = std::min(cal_upper_raw, precursor_mass_tolerance_upper_);
}
// else: cal_lower/cal_upper left at 0; writeback block must check !extreme_bias before use.
```

Updated writeback at cpp:770-780:
```cpp
if (cal.success && !cal.extreme_bias)
{
  Param fi_params = fi_params_original;
  fi_params.setValue("precursor:mass_tolerance_lower", cal.cal_lower);
  fi_params.setValue("precursor:mass_tolerance_upper", cal.cal_upper);
  fi_params.setValue("fragment:mass_tolerance", cal.fragment_tolerance);
  fragment_index_.setParameters(fi_params);
  fi_params_modified = true;

  // Refresh the algorithm's OWN member copies — OpenSearchModificationAnalysis (cpp:889-895)
  // reads them via computeModMatchTolerance_() and would otherwise see stale pre-calibration
  // values (latent bug today under the scalar tolerance; same class under the new schema).
  precursor_mass_tolerance_lower_ = cal.cal_lower;
  precursor_mass_tolerance_upper_ = cal.cal_upper;

  OPENMS_LOG_INFO << "[PDBS-FI] Calibration: shift=" << cal.precursor_shift
                  << " spread=" << cal.precursor_spread << " "
                  << precursor_mass_tolerance_unit_
                  << " → window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]" << std::endl;
}
else if (cal.success && cal.extreme_bias)
{
  OPENMS_LOG_WARN << "[PDBS-FI] Calibration: |shift|=" << std::abs(cal.precursor_shift)
                  << " >= spread=" << cal.precursor_spread << " "
                  << precursor_mass_tolerance_unit_ << " — calibration result discarded. "
                  << "The true signed window [" << (cal.precursor_shift - cal.precursor_spread)
                  << ", " << (cal.precursor_shift + cal.precursor_spread)
                  << "] lies entirely on one side of zero; the positive-magnitude schema cannot "
                  << "express it without loosening. Fix external calibration, or configure "
                  << "mass_tolerance_lower/_upper manually." << std::endl;
  // Preserve user-configured bounds. fi_params_modified stays false.
}
// cal.fragment_tolerance still applied separately when success — unchanged from pre-refactor.
```

**Why skip writeback on extreme bias.** The round-2 review worked the math: with `shift=+7, spread=2` the true calibrated window is `[+5, +9]`. Any positive-magnitude approximation either (a) floors `cal_lower` at 0 → effective window `[0, 9]`, which is **5 units wider** on the lower side than calibration justifies — a loosening, violating the "only tighten" invariant; or (b) encodes `[+5, +9]` exactly by storing `lower = -5` — breaks the positive-magnitude schema. Neither is acceptable. Discarding the calibration result and preserving user intent is the only honest choice, and matches the spec's "fix external calibration" advice.

**`!open_search` guard** at cpp:764 still works. `open_search` now comes from the new `isOpenSearchMode_()` delegating to the shared helper: `[20, 200]` Da → max 200 > 1 → open → calibration skipped. `[500, 1500]` ppm → max 1500 > 1000 → open → skipped.

**Fragment writeback unchanged.** `effective_fragment_tol` handling at cpp:773 is preserved — the new block replaces only the precursor-tolerance writeback.

### `OpenSearchModificationAnalysis` input (cpp:889-895) + new `computeModMatchTolerance_()` helper

`analyzeModifications` takes a scalar matching precision for Δmass → UniMod lookup. Under asymmetric bounds, the correct reduction is `min(lower, upper)` (tightest available matching precision). `OpenSearchModificationAnalysis.cpp:161-163` clamps every input at `MAX_MOD_MAPPING_TOL_ = 0.02 Da` and **discards ppm inputs entirely** (using the cap directly), so any realistic input is behavior-equivalent today — but `min` is the semantically honest choice and preserves correctness if the clamp is ever lifted.

Extract the reduction into a named private const helper so the test suite can observe it without spying on the `analyzeModifications` call (OpenMS class-tests have no spy/mock infrastructure):

```cpp
// PeptideSearchEngineFIAlgorithm.h (private section)
/// Scalar tolerance passed to OpenSearchModificationAnalysis under asymmetric bounds.
/// Uses the tighter of the two positive magnitudes — semantically correct for UniMod
/// matching precision. See spec §7 for the 0.02 Da internal clamp context.
double computeModMatchTolerance_() const
{
  return std::min(precursor_mass_tolerance_lower_, precursor_mass_tolerance_upper_);
}
```

Call site at cpp:889-895:
```cpp
mod_analyzer.analyzeModifications(
  peptide_ids,
  computeModMatchTolerance_(),
  precursor_mass_tolerance_unit_ == "ppm",
  false,
  "");
```

The helper is exposed to the test suite via a `PeptideSearchEngineFIAlgorithm_test` friend subclass in the test file (same pattern as `FragmentIndex_test` at `FragmentIndex_test.cpp:44`), so tests 8 and 10 in §9 can pin the reduction rule directly without running a full search or capturing log output.

## 8. Call-site changes: `PeptideDataBaseSearchFI.cpp`

The TOPP tool wrapper reads `precursor:mass_tolerance` and `mass_tolerance_unit` directly at `PeptideDataBaseSearchFI.cpp:168-169` and recomputes open-search mode locally at `:170-172`. Update both to read the new pair and apply the new detection rule. No other touches — the tool copies the full `Search:` param subtree to the algorithm and lets the algorithm handle the rest.

## 9. Tests

### Test infrastructure additions

OpenMS class-tests (`START_TEST` / `TEST_EQUAL`) have no mock/spy infrastructure. Observability for tests 7, 8, 9, 10 below uses named state, not intercepted calls. Tests 4 and 5 use an observable-proxy approach instead (running the same fixture twice with different settings and asserting PSM-set equivalence) — no new API required. The refactor introduces exactly four observability hooks, all minimal:

1. **`PeptideSearchEngineFIAlgorithm::computeModMatchTolerance_() const`** — private helper, returns `std::min(_lower_, _upper_)`. Already declared in §7 for code-structure reasons; also used as the test hook for the OSMA reduction rule. Exposed to tests via `using` in the test friend subclass.
2. **`CalibrationResult_::bool extreme_bias`** — already declared in §7. Test asserts on struct field instead of capturing `OPENMS_LOG_WARN`.
3. **`PeptideSearchEngineFIAlgorithm_test` friend subclass** (new, in `PeptideSearchEngineFIAlgorithm_test.cpp`) mirroring the `FragmentIndex_test` pattern at `FragmentIndex_test.cpp:44`. Exposes `precursor_mass_tolerance_lower_/_upper_`, `computeModMatchTolerance_`, `fragment_index_`, and `last_calibration_result_` (next hook).
4. **`PeptideSearchEngineFIAlgorithm::last_calibration_result_`** — private member, holds the most recent `CalibrationResult_` so tests can read `shift`, `spread`, `cal_lower`, `cal_upper`, `extreme_bias` after a search completes. Currently the calibration result is a local in the caller and is discarded; storing it is a one-line change. Also useful for debugging.

### Migration (mechanical)

- `src/tests/topp/PeptideDataBaseSearchFI_1.ini`: `mass_tolerance = 5` → `mass_tolerance_lower = 5`, `mass_tolerance_upper = 5`.
- `src/tests/topp/PeptideDataBaseSearchFI_2.ini`: `mass_tolerance = 10.0` → both bounds `10.0`.
- `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`: every `precursor:mass_tolerance` set site → positive-magnitude pair.
- `src/pyOpenMS/tests/test_modification_discovery.py`: **8 hits across 4 test functions** at lines 197-198, 300-301, 354-355, 434-435 (each function sets `mass_tolerance` + `mass_tolerance_unit`). **Not optional** — CI will fail on the first run otherwise. Note: pyOpenMS bindings are hand-maintained nanobind files (per `src/pyOpenMS/CLAUDE.md`), no codegen step — only the test file needs touching.

**Reference-file invariance (for symmetric configs).** `Math::ppmToMass<float>` is linear (`MathFunctions.h:403-407`): `ppmToMass<float>(10, m) = -ppmToMass<float>(-10, m)` in IEEE float (unary negation of a single product). The new code path computes `precursor_mass + {-Math::ppmToMass<float>(lower, m), +Math::ppmToMass<float>(upper, m)}` — bitwise identical to the legacy `precursor_mass ± Math::ppmToMass<float>(tol, m)` when `lower == upper == tol`. Any reference drift after migration is a real regression to root-cause, not a reference to bulk-update. The `static_cast<float>` chain on the `double` members and the explicit `<float>` template tag are load-bearing for this invariance — see §6.

### New tests (`FragmentIndex_test.cpp`)

1. **Symmetric default self-hit.** At default `[20, 20]` ppm, build the index from a small synthesized peptide set and confirm each peptide retrieves itself via `getPeptidesInMassWindow` within the expected range. Uses the existing `FragmentIndex_test` friend subclass pattern from `FragmentIndex_test.cpp:44-160`. Also confirm the pre-refactor code path on the same fixture yields the same peptide-index range — this is a bitwise-equivalence guard, not a stored-baseline comparison.
2. **Asymmetric window — simulate calibration offset.** Synthesize a spectrum whose `Precursor::getMZ()` is shifted by `setMZ(peptide.getMZ(1) * (1 + 8e-6))` (+8 ppm). At `[5, 5]` ppm, the peptide is NOT retrieved; at `[5, 15]` ppm it IS. Arithmetic check: at peptide mass 1500 Da, effective window `[-5, +15]` ppm = `[-0.0075, +0.0225]` Da; observed Δ = `+0.012` Da, which lies in `[5, 15]` but not in `[5, 5]` (`|Δ| > 0.0075`). Add this arithmetic as an inline comment so a future widening of the fixture doesn't silently break the boundary.
3. **Open-mode auto-detection boundaries** — call the static `FragmentIndex::isOpenSearchMode(lower, upper, unit_ppm)` directly (no `updateMembers_()` round-trip needed; the static is a pure function):
   - `isOpenSearchMode(500, 1500, true) == true` (flipped by upper)
   - `isOpenSearchMode(999, 999, true) == false`
   - `isOpenSearchMode(1000, 1000, true) == false` (strict `>`, exactly at threshold → closed)
   - `isOpenSearchMode(1000.0001, 1000, true) == true` (just over)
   - `isOpenSearchMode(0.9, 0.9, false) == false`, `isOpenSearchMode(1.0, 1.0, false) == false`, `isOpenSearchMode(1.1, 0.5, false) == true` (Da unit boundaries).
4. **Open-mode isotope-error force-to-zero — observable-proxy.** Build a fixture that depends on isotope_error iteration (e.g., target peptide whose precursor m/z is misassigned one `13C` up). Configure `isotope_error_min=-2, max=+2` AND open-mode-triggering tolerance. Compare PSM set against the same fixture with `isotope_error_min=0, max=0` + open mode. Equal PSM sets → iteration was collapsed. No internal-state inspection required.
5. **Asymmetric + isotope_error interaction — observable PSM-count check.** `[5, 15]` ppm closed mode with `isotope_error_range=[-1, +2]` on a fixture with peptides at `mass`, `mass + 1.00336 Da`, `mass + 2*1.00336 Da`, `mass - 1.00336 Da`. Each isotope-shifted peptide should be found at its corresponding shift; the `-1.00336 Da` peptide AND the `+2*1.00336 Da` peptide both fall inside the iteration. No white-box access needed — just candidate-set membership assertions.
6a. **Validation: negative magnitude** — `TEST_EXCEPTION(Exception::InvalidParameter, algo.setParameters(p_with_neg_lower))`. Fires from `Param::checkDefaults_` via `setMinFloat(0.0)` at the `setParameters` boundary, before the algorithm's `updateMembers_()` validator.
6b. **Validation: zero width** — `TEST_EXCEPTION(Exception::InvalidParameter, algo.setParameters(p_with_zero_both))`. Fires from `FragmentIndex::updateMembers_()` (where the explicit `lower + upper > 0` check lives), because `setMinFloat(0.0)` accepts zero.
6c. **Validation: NaN rejection** — `TEST_EXCEPTION(Exception::InvalidParameter, algo.setParameters(p_with_nan_upper))`. Fires from the `std::isfinite` guard in `FragmentIndex::updateMembers_()`.

### New tests (`PeptideSearchEngineFIAlgorithm_test.cpp`)

7. **Calibration preserves asymmetric bias (normal case).** Friend subclass exposes members. Configure `[20, 5]` ppm + `calibration:enabled=true`. Synthesize spectra (pattern: closed-search baseline around `PeptideSearchEngineFIAlgorithm_test.cpp:632+`) with a deliberate `+7 ppm` systematic shift on precursor m/z values. Run the search. Assert:
   - `algo.last_calibration_result_.success == true`
   - `algo.last_calibration_result_.extreme_bias == false`
   - `algo.last_calibration_result_.precursor_shift ≈ +7 ppm` (`TEST_REAL_SIMILAR` with loose tolerance)
   - `algo.last_calibration_result_.cal_lower < 20` and `cal_upper < 5` (tightened, not clobbered to symmetric)
   - `algo.precursor_mass_tolerance_lower_ == algo.last_calibration_result_.cal_lower` (member refresh fired)
   - `algo.fragment_index_.getParameters().getValue("precursor:mass_tolerance_lower") == cal_lower` (fragment-index side also updated)
8. **`computeModMatchTolerance_()` reads refreshed members.** Same setup as (7). After the search, assert `algo.computeModMatchTolerance_() == std::min(cal_lower, cal_upper)`. This is the double-bookkeeping regression guard — without the member-refresh lines in §7's calibration writeback, this would return `min(20, 5) == 5`, not the calibrated value. Pure helper call — no spy needed.
9. **Extreme bias path discards calibration.** Synthesize a fixture where `|shift| > spread` (e.g., all errors clustered near `+12 ppm` with very low variance). Assert:
   - `algo.last_calibration_result_.success == true`
   - `algo.last_calibration_result_.extreme_bias == true`
   - `algo.precursor_mass_tolerance_lower_/_upper_` equal the user-configured values (unchanged)
   - `algo.fragment_index_.getParameters()` also unchanged on the precursor side
   No log-stream capture — the bool flag is the test hook.
10. **Reduction rule pinned without calibration.** Disable calibration. Configure `[5, 50]` ppm → `computeModMatchTolerance_() == 5.0`. Configure `[50, 5]` ppm → `== 5.0`. Configure `[5, 50]` Da → `== 5.0`. Pure unit test on the helper via the friend subclass. Prevents silent drift if someone changes the reduction to `max` or midpoint.

### Regression watch

Run the full `FragmentIndex_test`, `PeptideSearchEngineFIAlgorithm_test`, the two TOPP `PeptideDataBaseSearchFI_*` integration tests, and `test_modification_discovery.py` after the refactor. Any diff in idXML reference outputs is a real regression to investigate, not a reference to bulk-accept.

## 10. Files touched

| File | Change |
|---|---|
| `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h` | rename `getPeptidesInPrecursorRange` → `getPeptidesInMassWindow` (const); add `computeMassWindow_` (const); new static `isOpenSearchMode(lower, upper, unit_ppm)`; update instance `isOpenSearchMode_()` to delegate; replace members with positive-magnitude pair; add `#include <algorithm>` |
| `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` | implementations of above; rewrite `searchDifferentPrecursorRanges`; `updateMembers_()` validation + WARN log on mode flip; param registrations mirror (drop `mass_tolerance`, `open_window_*`; add `mass_tolerance_lower/_upper`); delete dead "modification window" comment at `:993` |
| `src/openms/include/OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h` | replace scalar tolerance members with `_lower_` / `_upper_`; add private `computeModMatchTolerance_() const` helper; add `last_calibration_result_` member (stores the most recent `CalibrationResult_` for test observability); extend `CalibrationResult_` with `precursor_shift`, `precursor_spread`, `cal_lower`, `cal_upper`, `bool extreme_bias` (replace `precursor_tolerance`); `isOpenSearchMode_()` delegates to `FragmentIndex::isOpenSearchMode` static |
| `src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp` | param registrations (add pair, delete dead `open_window_*`, delete `mass_tolerance`); section description update; `updateMembers_()` reads the new pair; `runCalibrationPass_` signed-shift preservation with extreme-bias detection (cpp:1454-1626); calibration writeback skips on extreme bias, refreshes both `fragment_index_` params AND algo member copies in the normal case; `OpenSearchModificationAnalysis` input routes through `computeModMatchTolerance_()` |
| `src/topp/PeptideDataBaseSearchFI.cpp` | param reads (`mass_tolerance_lower/_upper`) + local open-mode detection (delegate to shared helper) |
| `src/tests/topp/PeptideDataBaseSearchFI_1.ini` / `_2.ini` | migrate param names (positive magnitudes) |
| `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp` | migrate + add `PeptideSearchEngineFIAlgorithm_test` friend subclass + new tests 7-10 (calibration normal case, helper refresh, extreme-bias path, reduction pinning) |
| `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` | add tests 1-6c above (symmetric self-hit, asymmetric fixture, boundary cases, observable-proxy iso tests, validation sub-cases incl. NaN) |
| `src/pyOpenMS/tests/test_modification_discovery.py` | **migrate 8 hits across 4 test functions** at lines 197-198, 300-301, 354-355, 434-435 — **mandatory** (otherwise first CI run red). pyOpenMS bindings are hand-maintained; no codegen step |
| `CHANGELOG` | `BREAKING:` entry under `TOPP tools: > Changes: > PeptideDataBaseSearchFI:` with `(#9108)` suffix, noting param renames + pyOpenMS user impact |

## 11. Deferred / follow-ups

- **Wide-window acquisition (WWA) implementation.** NOT in this PR — no flag, no param, no stub. The `getPeptidesInMassWindow` candidate-lookup *signature* is stable and will accept WWA inputs without further change. What the WWA PR will need to add:
  1. **`precursor:wide_window` bool param + member.** Declared for the first time in the WWA PR, where the name can be validated against the actual implementation rather than committed speculatively now.
  2. **Per-spectrum Th→Da conversion** in `PeptideSearchEngineFIAlgorithm`. `Precursor::getIsolationWindowLowerOffset()` / `UpperOffset()` return m/z (Thomson) — the caller must multiply by the assigned charge to get Da offsets before handing them to `getPeptidesInMassWindow(mass, {lo_Da, hi_Da})`.
  3. **Charge-loop branch** in `FragmentIndex::querySpectrum()` at `FragmentIndex.cpp:1027-1054`. Today the code short-circuits to the single reported charge when `precursor.getCharge() != 0` (`:1030-1033`). WWA mode must override that short-circuit and iterate the plausible-charge range, because the wide isolation window contains co-fragmenting species at **other** charges.
  4. **Isotope-error force-to-zero** under wide-window, same as open mode (the isotope offset is redundant with a wide isolation window).
  5. **Chimera gate at `FragmentIndex.cpp:1021`** (`precursor.size() != 1`) is **preserved** by this refactor and intentionally remains compatible with WWA — WWA DDA spectra still list one assigned precursor with a wide isolation window.
- **Dedicated mod-matching tolerance.** `modification:matching_tolerance` param for `OpenSearchModificationAnalysis`, decoupled from the (wide) search window and from the internal `MAX_MOD_MAPPING_TOL_ = 0.02 Da` clamp. Improves UniMod matching precision when the clamp is lifted. Separate design.
- **Lifting `MAX_MOD_MAPPING_TOL_` clamp.** Currently hidden inside `OpenSearchModificationAnalysis.cpp:161-163`; the ppm path discards its input entirely. Revisit as part of the dedicated mod-matching tolerance work.
- **Completing `const`-correctness on `FragmentIndex`.** `queryPeaks`, `query`, `searchDifferentPrecursorRanges`, `querySpectrum` are all non-`const` today and none of them actually mutate the index. Out of scope here; worth a follow-up janitor pass.

## 12. Risk

- **Medium** — driven by the calibration rewrite. The pure-refactor portion (asymmetric-param rename + `FragmentIndex` interface refactor) is low risk; the calibration bias-preservation work (§7) is new functionality with new math, a new struct shape, and new edge-case logic (extreme-bias detection and writeback skip). Tests 7, 9, 10 in §9 exercise it directly.
- **pyOpenMS breakage surface.** `PeptideSearchEngineFIAlgorithm` is fully bound via `bind_misc.cpp:3322-3406`. Python users hitting the old `precursor:mass_tolerance` param name will get a `checkDefaults_` error the first time they call `setParameters()`. Mitigation: mandatory `test_modification_discovery.py` migration (§10) + CHANGELOG `BREAKING:` entry naming the rename explicitly.
- **Test reference drift** — mitigated by the bitwise-equivalence argument in §9 plus the explicit "any diff is a real regression" rule. The `(float)` cast chain + explicit `<float>` template tag are load-bearing; see §6 and §9.
- **Calibration double-bookkeeping.** The fix writes to two places (`fragment_index_` params AND algo-level members). Missing either side reintroduces the latent stale-tolerance bug that already exists today under the scalar schema. Tests 7 and 8 in §9 pin both sides independently via the friend subclass.
