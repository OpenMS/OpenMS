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

- **G1.** Users can specify an asymmetric precursor window for closed search (e.g., `mass_tolerance_lower=5, mass_tolerance_upper=15` ppm → effective window `[−5, +15]` ppm to compensate a known calibration bias). Calibration pass (§7) preserves this bias instead of clobbering it to symmetric.
- **G2.** Closed and open search use one unified asymmetric window mechanism. No duplicate knobs.
- **G3.** The `FragmentIndex` candidate-lookup API has a single, honest, absolute-bounds signature. The current additive `prec_tol + window` semantics — which hide the tolerance inside the function and only accept an *additional* offset from the caller — are removed.
- **G4.** `wide_window` flag declared (default `false`, implementation deferred) so the per-spectrum WWA path can wire into the same refactored interface without further churn.
- **G5.** No regression in modification discovery for today's open-search flow.

## 3. Non-goals

- **NG1.** Implementing wide-window acquisition / per-spectrum isolation-window candidate selection. The flag is declared; the implementation is a separate design.
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
- Users rarely need truly asymmetric windows in practice; magnitudes read naturally as "how much slack in each direction".
- Sign flips happen exactly once, inside `computeMassWindow_()`, where the convention is centralized.

## 5. Parameter schema

### Added

| Param | Type | Default | Notes |
|---|---|---|---|
| `precursor:mass_tolerance_lower` | `double` | `20.0` | **Positive magnitude.** Applied internally as `−lower` to form the lower bound. Must be ≥ 0. Unit controlled by `precursor:mass_tolerance_unit`. Default widened from the original `10.0` to reduce first-run surprises on uncalibrated / non-Orbitrap instruments (Bruker timsTOF typically needs 15–20 ppm). The internal calibration pass (§7) tightens it on real data. |
| `precursor:mass_tolerance_upper` | `double` | `20.0` | **Positive magnitude.** Applied internally as `+upper`. Must be ≥ 0. Same unit. |
| `precursor:wide_window` | `string` (enum `"true"`/`"false"`) | `"false"` | **Declared but not implemented in this PR.** Registered only in `PeptideSearchEngineFIAlgorithm` (not `FragmentIndex` — see §7), since the flag requires per-spectrum metadata the index never touches. When `true`, the algorithm throws `Exception::NotImplemented` at `updateMembers_()`-time. Future: the WWA path reads per-spectrum isolation-window offsets and bypasses the member bounds. |

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

- Both `mass_tolerance_lower` and `mass_tolerance_upper` must be `>= 0` (positive magnitudes). Nonzero total width (`lower + upper > 0`) required.
- Unit `"Da"` accepts any magnitude; unit `"ppm"` is unbounded but auto-trips open-mode detection below.

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
// NB: no wide_window_ member here — lives in PeptideSearchEngineFIAlgorithm (§7).
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

// Validation
if (precursor_mass_tolerance_lower_ < 0.0 || precursor_mass_tolerance_upper_ < 0.0)
{
  throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    "precursor:mass_tolerance_lower and mass_tolerance_upper must be >= 0 (positive magnitudes)");
}
if (precursor_mass_tolerance_lower_ + precursor_mass_tolerance_upper_ <= 0.0)
{
  throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    "precursor window has zero width");
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

`wide_window` is intentionally **not** read here — it's a `PeptideSearchEngineFIAlgorithm`-layer concern (§7).

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

1. **Dead parameter registration** (`PeptideSearchEngineFIAlgorithm.cpp:181-182`): registers `precursor:open_window_lower`/`_upper` but never reads them in `updateMembers_()`. Pure pass-through to `FragmentIndex`. **Delete these two lines.**
2. **Duplicate `isOpenSearchMode_()`** (`PeptideSearchEngineFIAlgorithm.h:478-483`): a second copy of the auto-detection logic lives as an inline method in the header, alongside `FragmentIndex`'s. Consolidate via the new static helper on `FragmentIndex` (below).
3. **`OpenSearchModificationAnalysis` tolerance input** (`PeptideSearchEngineFIAlgorithm.cpp:889-895`): currently passes `precursor_mass_tolerance_` (the scalar) and the unit. It uses this tolerance to match observed Δmass to UniMod entries.
4. **Calibration pass discards signed bias** (`PeptideSearchEngineFIAlgorithm.cpp:1580`): `runCalibrationPass_` computes signed `prec_err` at `:1534` then throws the sign away via `std::abs` at `:1580`. Also updates **only** the `fragment_index_` params, leaving the algo-level `precursor_mass_tolerance_*_` stale — a latent bug today (since `OpenSearchModificationAnalysis` reads the stale members), which would silently survive the refactor. Both fixed here.

### Parameter registration changes (PeptideSearchEngineFIAlgorithm.cpp:65-72, 181-182)

Replace
```cpp
defaults_.setValue("precursor:mass_tolerance", 10.0, "+/- tolerance for precursor mass.");
```
with
```cpp
defaults_.setValue("precursor:mass_tolerance_lower", 20.0,
                   "Lower-side precursor-mass slack (positive magnitude, applied as -lower). "
                   "Note: when strongly asymmetric, also review precursor:isotope_error_min.");
defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
defaults_.setValue("precursor:mass_tolerance_upper", 20.0,
                   "Upper-side precursor-mass slack (positive magnitude, applied as +upper).");
defaults_.setMinFloat("precursor:mass_tolerance_upper", 0.0);
defaults_.setValue("precursor:wide_window", "false",
                   "Declared for forward compatibility — per-spectrum isolation-window search. "
                   "NOT YET IMPLEMENTED; setting to true throws Exception::NotImplemented.");
defaults_.setValidStrings("precursor:wide_window", {"true", "false"});
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
wide_window_                    = param_.getValue("precursor:wide_window").toString() == "true";

if (wide_window_)
{
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
}
```

`wide_window_` is a new private `bool` member in `PeptideSearchEngineFIAlgorithm.h`, declared next to the existing precursor-tolerance members. `FragmentIndex` never sees the flag.

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

### Calibration pass changes — preserve signed bias (cpp:760-782, 1454-1626)

**Blocker fix from review.** The current `runCalibrationPass_` computes signed `prec_err` at cpp:1534, then immediately `std::abs()`-es it at cpp:1580 and returns a single-scalar spread. Writing that scalar back as `[−spread, +spread]` discards any real calibration bias — which is precisely the use case that motivates asymmetric closed search. A user setting `[15, 5]` ppm to compensate a −5 ppm bias would have the window clobbered to `[3, 3]` ppm and lose every target.

**Fix:** extend `CalibrationResult_` (`PeptideSearchEngineFIAlgorithm.h:444-450`) with a signed shift, compute `median(signed_err)` + MAD-around-the-median, and write asymmetric bounds back.

Extended struct:
```cpp
struct CalibrationResult_
{
  double precursor_shift{0};    ///< median of signed precursor errors (bias, same unit as configured)
  double precursor_spread{0};   ///< half-width around the shift (median + 3*MAD of |err - shift|)
  double fragment_tolerance{0};
  double fragment_shift{0};
  bool   success{false};
};
```
(`precursor_tolerance` is replaced by the more honest `precursor_shift` + `precursor_spread` pair.)

Updated `runCalibrationPass_` (replacing cpp:1575-1601):
```cpp
// Signed median gives the calibration bias; MAD around it gives the spread.
std::vector<double> signed_errors = precursor_errors; // already signed from cpp:1534-1536
const double prec_shift = Math::median(signed_errors.begin(), signed_errors.end());

std::vector<double> residuals;
residuals.reserve(signed_errors.size());
for (double e : signed_errors) { residuals.push_back(std::abs(e - prec_shift)); }
const double prec_mad = Math::median(residuals.begin(), residuals.end());

result.precursor_shift  = prec_shift;
result.precursor_spread = std::max(min_tolerance, 3.0 * prec_mad);
```
The "don't widen beyond configured" cap at cpp:1605-1608 generalizes to asymmetric:
```cpp
// Only tighten. Cap each side independently against the user-configured bound.
const double cal_lower_raw = std::max(0.0, result.precursor_spread - result.precursor_shift);
const double cal_upper_raw = std::max(0.0, result.precursor_spread + result.precursor_shift);
result.cal_lower = std::min(cal_lower_raw, precursor_mass_tolerance_lower_);
result.cal_upper = std::min(cal_upper_raw, precursor_mass_tolerance_upper_);
```
(`cal_lower` and `cal_upper` are added to `CalibrationResult_` alongside `shift`/`spread`.)

Updated writeback at cpp:770-780:
```cpp
if (cal.success)
{
  Param fi_params = fi_params_original;
  fi_params.setValue("precursor:mass_tolerance_lower", cal.cal_lower);
  fi_params.setValue("precursor:mass_tolerance_upper", cal.cal_upper);
  fi_params.setValue("fragment:mass_tolerance", cal.fragment_tolerance);
  fragment_index_.setParameters(fi_params);
  fi_params_modified = true;

  // **Fix for the double-bookkeeping blocker**: refresh the algorithm's OWN member copies too.
  // Otherwise OpenSearchModificationAnalysis (cpp:889-895) reads stale pre-calibration values.
  precursor_mass_tolerance_lower_ = cal.cal_lower;
  precursor_mass_tolerance_upper_ = cal.cal_upper;
  effective_fragment_tol           = cal.fragment_tolerance;

  OPENMS_LOG_INFO << "[PDBS-FI] Calibration: shift=" << cal.precursor_shift
                  << " spread=" << cal.precursor_spread << " "
                  << precursor_mass_tolerance_unit_
                  << " → window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]" << std::endl;
  if (cal.precursor_shift > cal.precursor_spread
      || -cal.precursor_shift > cal.precursor_spread)
  {
    OPENMS_LOG_WARN << "[PDBS-FI] Calibration: shift exceeds spread (extreme instrument bias). "
                    << "Window floored at one edge; consider fixing calibration externally "
                    << "before running this tool." << std::endl;
  }
}
```

**Why the two member-refresh lines matter.** Today, `OpenSearchModificationAnalysis` at cpp:889-895 reads `precursor_mass_tolerance_` from the algorithm's members — which, after calibration, is **stale** (FragmentIndex got updated, the algo didn't). Under the new asymmetric flow, the same bug would fire via `precursor_mass_tolerance_lower_/_upper_`. This refactor fixes it alongside the calibration rewrite rather than leaving it as pre-existing silent breakage.

**Positive-magnitude edge case.** If `|shift| > spread`, the signed window `[shift−spread, shift+spread]` lies entirely on one side of zero — which the positive-magnitude schema can't express, since `[-lower, +upper]` always straddles zero. In that case we floor the opposite side at 0 and emit a `WARN`. This is the exact point where the user should fix external calibration; the search-engine should not silently mask it. The same convention is used by `MSFraggerAdapter`.

**`!open_search` guard** at cpp:764. `open_search` is now computed from the new `isOpenSearchMode_()` (which in turn delegates to the shared helper). `[−100, +200]` Da → `max(100,200) > 1` → open → calibration skipped. `[−500, +1500]` ppm → `max(500,1500) > 1000` → open → skipped. Correct.

### `OpenSearchModificationAnalysis` input (cpp:889-895)

**Corrected from the original spec draft after review.** The scalar passed to `analyzeModifications` is a **matching precision** for Δmass → UniMod lookup — tighter = better, not wider = safer. The original draft used `max(|lower|, upper)`, which would *loosen* matching for asymmetric users. Fix: use `min(lower, upper)`.

However, reading `OpenSearchModificationAnalysis.cpp:155-163`:
```cpp
const double effective_tol = precursor_mass_tolerance_unit_ppm
  ? MAX_MOD_MAPPING_TOL_                                    // 0.02 Da, ignores ppm input entirely
  : std::min(precursor_mass_tolerance, MAX_MOD_MAPPING_TOL_); // Da: hard cap at 0.02 Da
```
reveals that the scalar is **clamped to `MAX_MOD_MAPPING_TOL_ = 0.02 Da` unconditionally**. In ppm mode the input value is **discarded entirely**. So in practice `min`, `max`, midpoint, or the original scalar all produce bitwise-identical output for every realistic configuration (any ppm value → 0.02 Da; any Da value ≥ 0.02 → 0.02 Da).

Pick `min(lower, upper)` anyway, because:
1. It is semantically correct (matching precision).
2. It is strictly ≤ today's symmetric scalar, so no regression is possible even if the clamp is later lifted.
3. The comment documents the invariant explicitly so a future maintainer "fixing" the clamp doesn't accidentally widen matching.

```cpp
// OpenSearchModificationAnalysis matches observed Δmass to UniMod entries using a
// scalar matching precision (tighter = better). Under asymmetric bounds we feed the
// tighter of the two magnitudes — strictly ≤ today's symmetric value, never worse.
//
// Note: internally, analyzeModifications clamps this at MAX_MOD_MAPPING_TOL_ = 0.02 Da
// (and DISCARDS the value entirely in ppm mode, using the cap directly). So in
// practice any realistic input produces the same output. This is correct behavior
// today but a latent regression pit if the clamp is ever lifted without revisiting
// this reduction rule. See spec §11 for the follow-up "dedicated mod-matching
// tolerance" item.
const double mod_match_tol = std::min(precursor_mass_tolerance_lower_,
                                      precursor_mass_tolerance_upper_);
mod_analyzer.analyzeModifications(
  peptide_ids,
  mod_match_tol,
  precursor_mass_tolerance_unit_ == "ppm",
  false,
  "");
```

No change to `OpenSearchModificationAnalysis` itself. Testing: §9 adds an explicit regression test pinning the input value that reaches `analyzeModifications`, so drift in this reduction rule is detectable even while the internal clamp absorbs the effect.

## 8. Call-site changes: `PeptideDataBaseSearchFI.cpp`

The TOPP tool wrapper reads `precursor:mass_tolerance` and `mass_tolerance_unit` directly at `PeptideDataBaseSearchFI.cpp:168-169` and recomputes open-search mode locally at `:170-172`. Update both to read the new pair and apply the new detection rule. No other touches — the tool copies the full `Search:` param subtree to the algorithm and lets the algorithm handle the rest.

## 9. Tests

### Migration (mechanical)

- `src/tests/topp/PeptideDataBaseSearchFI_1.ini`: `mass_tolerance = 5` → `mass_tolerance_lower = 5`, `mass_tolerance_upper = 5` (positive magnitudes).
- `src/tests/topp/PeptideDataBaseSearchFI_2.ini`: `mass_tolerance = 10.0` → `mass_tolerance_lower = 10.0`, `mass_tolerance_upper = 10.0`.
- `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp`: every site that sets `precursor:mass_tolerance` → replace with the positive-magnitude pair.
- `src/pyOpenMS/tests/test_modification_discovery.py`: `precursor:mass_tolerance` at lines 197-198, 300-301, 354-355, 434-435 (4 call sites) → migrate in lockstep. **Not optional** — this test file will fail on the first CI run after the C++ change otherwise.

**Reference-file invariance claim.** idXML reference outputs for the migrated symmetric INIs should be **unchanged**. Proof sketch: `Math::ppmToMass` is linear (`MathFunctions.h:403-407`), so `Math::ppmToMass<float>(10, m)` and `-Math::ppmToMass<float>(10, m)` are bitwise additive inverses in IEEE float. The new code path `precursor_mass + {-Math::ppmToMass<float>(lower, m), +Math::ppmToMass<float>(upper, m)}` produces the same float values that the old path computed as `precursor_mass ± Math::ppmToMass<float>(tol, m)`, modulo the `static_cast<float>` chain on the member inputs (which matters — see "load-bearing" note below). If any reference file diffs after migration, it's a real regression to root-cause, not a reference to bulk-update.

**`(float)` cast chain is load-bearing.** The equivalence above relies on `computeMassWindow_` casting the `double` members to `float` *before* calling `Math::ppmToMass<float>(...)`, exactly matching how the legacy code computed `prec_tol` as `float`. If someone later widens the helper to `double`, the final `float`-narrowing of the final mass comparison differs and reference files drift. Explicit `<float>` template tag + `static_cast<float>` on the member input guards against this.

### New tests (`FragmentIndex_test.cpp`)

1. **Symmetric default equivalence.** Build index at default `[20, 20]` ppm, search a fixture spectrum, confirm the PSM set is identical to the pre-refactor baseline. Guards the refactor via bitwise `lower_bound` invariance.
2. **Asymmetric closed window — simulate calibration offset.** Load a fixture spectrum, **apply `Precursor::setMZ(exp_mz * (1 + 8e-6))` to shift the *spectrum precursor* by +8 ppm** (the DB is correct, the instrument is mis-calibrated — the realistic scenario). With `[5, 5]` ppm bounds the target peptide is NOT identified; with `[5, 15]` ppm it IS. Template pattern: existing `isotope_error` section (`FragmentIndex_test.cpp:453-455`).
3. **Open-mode auto-detection boundaries.**
   - `[500, 1500]` ppm → open mode (flipped by upper side, strict `>` at 1000 threshold)
   - `[999, 999]` ppm → closed mode
   - `[1000, 1000]` ppm → closed mode (strict `>`, exactly at threshold)
   - `[1000.0001, 1000]` ppm → open mode (just over)
4. **Open-mode isotope-error force-to-zero.** Confirm iteration collapses to `[0, 0]` when open mode trips, even with user-configured `isotope_error_min=−2, max=+2`.
5. **Asymmetric + isotope-error interaction.** `[5, 15]` ppm closed mode with `isotope_error_range = [-1, +2]`. Verify the window travels correctly with each `+k·1.00336 Da` shift and candidate sets match the expected union.
6. **Validation.** `mass_tolerance_lower = -5` throws `Exception::InvalidParameter`. `[0, 0]` throws `Exception::InvalidParameter` (zero width). Use `TEST_EXCEPTION(Exception::InvalidParameter, ...)` macro.
7. **`wide_window = "true"` throws `Exception::NotImplemented`.** Declared in `PeptideSearchEngineFIAlgorithm`, test lives in `PeptideSearchEngineFIAlgorithm_test.cpp`. Use `TEST_EXCEPTION(Exception::NotImplemented, ...)`.

### New tests (`PeptideSearchEngineFIAlgorithm_test.cpp`)

8. **Calibration preserves asymmetric bias.** Configure `[20, 5]` ppm with a fixture run whose signed `prec_err` has `median ≈ +7 ppm` (synthesized or captured). Run the algorithm with `calibration:enabled=true`. After the search, read back `fragment_index_.getParameters()` and the algo's own `precursor_mass_tolerance_lower_/_upper_` members. Assert both sides reflect the calibrated `[cal_lower, cal_upper]` (not the original `[20, 5]`), and that neither has been clobbered to the symmetric scalar of previous behavior. Capture the `OPENMS_LOG_INFO` "Calibration: shift=... spread=..." line to verify the signed shift is reported.
9. **Algorithm-level member refresh after calibration.** Same setup as (8). Spy on the `OpenSearchModificationAnalysis::analyzeModifications` call at cpp:889-895 and assert the `precursor_mass_tolerance` argument equals `min(cal.cal_lower, cal.cal_upper)` — NOT the user-configured `min(20, 5) = 5`. This is the double-bookkeeping fix regression guard.
10. **Calibration: extreme bias warning path.** Construct a fixture where `|shift| > spread`. Assert `OPENMS_LOG_WARN` fires and one side of the window is floored at 0.
11. **OpenSearchModificationAnalysis reduction under asymmetry (no calibration).** `[5, 50]` ppm and `[50, 5]` ppm → both feed `5` ppm into `analyzeModifications`. `[5, 50]` Da → feeds 5 Da. Pins the `min()` reduction rule even though the internal 0.02 Da clamp absorbs the effect. Prevents silent drift if someone changes the reduction to `max` or midpoint.

### Regression watch

Run the full `FragmentIndex` and `PeptideSearchEngineFIAlgorithm` class tests + the two TOPP `PeptideDataBaseSearchFI_*` integration tests + `test_modification_discovery.py` after the refactor. Any diff in reference outputs is a real regression, not a reference to update.

## 10. Files touched

| File | Change |
|---|---|
| `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h` | rename `getPeptidesInPrecursorRange` → `getPeptidesInMassWindow` (const); add `computeMassWindow_` (const); new static `isOpenSearchMode(lower, upper, unit_ppm)`; update instance `isOpenSearchMode_()` to delegate; replace members with positive-magnitude pair; add `#include <algorithm>` |
| `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` | implementations of above; rewrite `searchDifferentPrecursorRanges`; `updateMembers_()` validation + WARN log on mode flip; param registrations mirror (drop `mass_tolerance`, `open_window_*`; add `mass_tolerance_lower/_upper`); delete dead "modification window" comment at `:993` |
| `src/openms/include/OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h` | replace scalar tolerance members with `_lower_` / `_upper_`; add `wide_window_` bool member; extend `CalibrationResult_` with `precursor_shift`, `precursor_spread`, `cal_lower`, `cal_upper` (replace `precursor_tolerance`); `isOpenSearchMode_()` delegates to `FragmentIndex::isOpenSearchMode` static |
| `src/openms/source/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.cpp` | param registrations (add pair + `wide_window`, delete dead `open_window_*`, delete `mass_tolerance`); section description update; `updateMembers_()` with `wide_window` read + throw; `runCalibrationPass_` signed-shift preservation (lines 1520-1626); calibration writeback refreshes both `fragment_index_` params AND algo member copies; `OpenSearchModificationAnalysis` input uses `min(lower, upper)` + comment block |
| `src/topp/PeptideDataBaseSearchFI.cpp` | param reads (`mass_tolerance_lower/_upper`) + local open-mode detection (delegate to shared helper) |
| `src/tests/topp/PeptideDataBaseSearchFI_1.ini` / `_2.ini` | migrate param names (positive magnitudes) |
| `src/tests/class_tests/openms/source/PeptideSearchEngineFIAlgorithm_test.cpp` | migrate + new tests (8)-(11): calibration preservation, member refresh, extreme-bias warn, mod-analyzer reduction |
| `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` | add tests (1)-(7) above |
| `src/pyOpenMS/tests/test_modification_discovery.py` | **migrate param names** at lines 197-198, 300-301, 354-355, 434-435 — **mandatory** (otherwise first CI run red) |
| `CHANGELOG` | `BREAKING:` entry under `TOPP tools: > Changes: > PeptideDataBaseSearchFI:` with `(#9108)` suffix, noting param renames + pyOpenMS user impact |

## 11. Deferred / follow-ups

- **Wide-window implementation.** The `getPeptidesInMassWindow` candidate-lookup *signature* is stable. What the WWA PR will still need to add:
  1. **Per-spectrum Th→Da conversion** in `PeptideSearchEngineFIAlgorithm`. `Precursor::getIsolationWindowLowerOffset()` / `UpperOffset()` return m/z (Thomson) — the caller must multiply by the assigned charge to get Da offsets before handing them to `getPeptidesInMassWindow(mass, {lo_Da, hi_Da})`.
  2. **Charge-loop branch** in `FragmentIndex::querySpectrum()` at `FragmentIndex.cpp:1027-1054`. Today the code short-circuits to the single reported charge when `precursor.getCharge() != 0` (`:1030-1033`). WWA mode must override that short-circuit and iterate the plausible-charge range, because the wide isolation window contains co-fragmenting species at **other** charges.
  3. **Isotope-error force-to-zero** under wide-window, same as open mode (the isotope offset is redundant with a wide isolation window).
  4. **Chimera gate at `FragmentIndex.cpp:1021`** (`precursor.size() != 1`) is **preserved** by this refactor and intentionally remains compatible with WWA — WWA DDA spectra still list one assigned precursor with a wide isolation window.
- **Dedicated mod-matching tolerance.** `modification:matching_tolerance` param for `OpenSearchModificationAnalysis`, decoupled from the (wide) search window and from the internal `MAX_MOD_MAPPING_TOL_ = 0.02 Da` clamp. Improves UniMod matching precision when the clamp is lifted. Separate design.
- **Lifting `MAX_MOD_MAPPING_TOL_` clamp.** Currently hidden inside `OpenSearchModificationAnalysis.cpp:161-163`; the ppm path discards its input entirely. Revisit as part of the dedicated mod-matching tolerance work.
- **Completing `const`-correctness on `FragmentIndex`.** `queryPeaks`, `query`, `searchDifferentPrecursorRanges`, `querySpectrum` are all non-`const` today and none of them actually mutate the index. Out of scope here; worth a follow-up janitor pass.

## 12. Risk

- **Low-to-medium** — changes are confined to one search engine with no released users, BUT the calibration-bias preservation work (§7) genuinely extends `CalibrationResult_` and touches `runCalibrationPass_` internals. That's new functionality riding alongside the refactor, not a pure rename.
- **pyOpenMS breakage surface.** `PeptideSearchEngineFIAlgorithm` is fully bound via `bind_misc.cpp:3322-3406`. Python users hitting the old `precursor:mass_tolerance` param name will get silent "unknown key" behavior until `checkDefaults_()` fires. The mitigation is the mandatory `test_modification_discovery.py` migration (§10) + a CHANGELOG `BREAKING:` entry that names the rename explicitly.
- **Test reference drift** — mitigated by the bitwise-equivalence argument in §9 plus the explicit "any diff is a real regression, not a reference to bulk-update" rule. The `(float)` cast chain is load-bearing; see §9.
- **Calibration-refresh blast radius.** The fix writes back to two places (`fragment_index_` params AND algo-level members). Missing either side reintroduces the `OpenSearchModificationAnalysis` stale-tolerance bug that exists today. Tests 8 and 9 in §9 pin both sides independently.
