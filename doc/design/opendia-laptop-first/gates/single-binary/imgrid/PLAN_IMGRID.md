# PLAN_IMGRID (frozen before the 5-seed runs; precheck seed 0 on D was run, see PRECHECK_IMGRID.md)

## Change (on top of 6dd893b954)
- `FragmentGridCollector::enableIonMobility(tol)`: with a mobility calibration, the 12-fragment grid sums only points
  with |1/K0 − calibrated precursor 1/K0| ≤ the calibrated im_tolerance. Precursors without a library 1/K0, and
  spectra without VSSC mobility, are pooled as before.
- `InterferenceResolver::compute(..., im_link_width = 2 × im_tolerance)`: precursors whose calibrated 1/K0 differ by
  more than this do not link. With an unknown 1/K0 the link uses m/z only.
- Both are on by default exactly when the calibration has a mobility model. No new constant: the width is the run's
  own calibrated tolerance. Astral (no mobility) is unchanged.

## Runs
- BASE: 6dd893b954 (`OpenDIAWorkflow.rtcentre`). NEW: the commit with this change (`OpenDIAWorkflow.imgrid`).
- D: HeLa 2737.d, seeds 0-4, BASE runs from PLAN_RTCENTRE reused (same binary). NEW 5 seeds.
- H (held out): HeLa 22474.d, seeds 0-4, BASE from PLAN_RTCENTRE reused. NEW evaluated only after the D/A verdict is
  committed.
- A: Astral raw01 seed 0 NEW must be byte-identical to BASE (PLAN_RTCENTRE A_rtcentre).

## Metrics
N (own q ≤ 0.01) and agreement (share of N in DIA-NN's q ≤ 0.01 list), as in PLAN_RTCENTRE (`eval_rtc.py`).

## Pass rules (all must hold)
1. D: median N(NEW) ≥ 1.03 × median N(BASE), and min N(NEW) > max N(BASE).
2. D: median agreement(NEW) ≥ median agreement(BASE) − 0.005.
3. A: byte-identical.
4. H (after the verdict): median N(NEW) ≥ median N(BASE), and median agreement(NEW) ≥ median agreement(BASE) − 0.005.

On failure: revert, no retuning.
