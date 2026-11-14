# Calibration decoy sensitivity: diagnosis and fix (OpenMS 0e4a234762)

## Cause
- Stage 3 of OpenDIACalibrate kept anchors by a paired target-decoy competition on the co-elution score alone.
- With the old decoys (copied assays), decoys rarely reached high co-elution scores:
  - PAIRLIB raw01: 0 decoys among the top 100;
  - 607 anchors at q ≤ 0.01.
- With predicted decoys, 3–33 decoys reach the top 100, and stage 3 accepts 0 anchors (`calib/armC_raw01`, reproduced
  locally).
- The top decoys are mostly chance co-elution, not signal borrowed from their own target. Excluding decoys that share
  ≥ 1 or ≥ 2 fragments with their target does not reliably restore anchors (`overlap_sim.py`: 413, 0, 401).
- What separates targets from these decoys is RT. A chance peak lies uniformly over the stage-3 window, a target near
  the map. Combining the score's p-value among decoys with an RT likelihood ratio restores 372 anchors on arm C raw01
  offline (`rt_term_sim.py`). It also raises the old library to 750/806 anchors (from 607/661).

## Fix
Stage 3 ranks by −log p(score | panel decoys) + log N(residual; 0, stage-1 σ) / U(window). Both terms are measured on
the panel; there are no new parameters.

## Recalibration, all libraries (`calib/new_*`, `map_accuracy.txt`)
| run | before | after: anchors | map vs DIA-NN, abs. error p50 / p95 (s) | DIA-NN IDs within own 4σ |
|---|---|---|---|---|
| old library raw01 / raw02 | 607 / 661 | 753 / 807 | 17.4 / 85.4, 17.2 / 85.9 (before: 16.9 / 80.8, 16.9 / 85.2) | 0.943 / 0.938 (before: 0.957 / 0.946) |
| A raw01 / raw02 | 225 / 205 | 524 / 466 | 17.5 / 83.9, 17.3 / 85.1 | 0.932 / 0.903 |
| B raw01 / raw02 | fail / 122 | **fail** / 305 | –, 18.3 / 69.1 | – / 0.915 |
| C raw01 / raw02 | fail / fail | 374 / 413 | 17.6 / 69.9, 17.2 / 69.2 | 0.904 / 0.929 |
| D raw01 / raw02 | 515 / 371 | 455 / 662 | 16.4 / 63.6, 16.9 / 69.6 | 0.963 / 0.954 |

- **Fixed:** the stage-3 failures. Map accuracy is unchanged.
- **Not fixed:** arm B raw01 still fails in stage 1. That is a whole-run competition with no RT information, and its
  failure has a different cause.
- **Cost:** anchors chosen with an RT preference give a tighter map σ (13.5–18.9 s, against 19.8–22.7 s). As a result,
  4σ covers 90–96 % of DIA-NN IDs, against 95 %. The σ estimate is biased low by the selection. A fix (for example a σ
  estimated from the score-only anchors of the same window, or a truncated-Gaussian correction) is not done.
- The arm results in RESULT_DECOYGEN used the old PAIRLIB calibration (amendment 1) and are unaffected.

## Follow-up: selection-free RT spread (OpenMS commit after 0e4a234762, see git log)
- **Problem:** selecting anchors with an RT preference biases their spread low.
  - DIA-NN's robust spread around the same maps is 24.4–27.1 s.
  - The written σ was 13.5–18.9 s.
  - So 4σ covered only 90–96 % of DIA-NN IDs.
- **Fix:** stage 4 fits all stage-3 target winners' residuals to the final map as a Gaussian (true peaks) plus
  background, the background shape taken from the decoy winners' residuals over the same window. The tolerance is
  4 × that width (`rt.spread_s`). It needs no labels and no threshold, and falls back to the map σ without enough
  signal.
  - Offline check (`mixture_sigma.py`): the fit recovered DIA-NN's spread on every run, including the old library.

| run | map σ (anchors) | spread_s | tolerance | DIA-NN 1 % IDs within tolerance |
|---|---|---|---|---|
| old code, old library raw01 / raw02 | 22.71 / 19.84 | – | 90.8 / 79.4 s | 0.957 / 0.946 |
| old library raw01 / raw02 | 18.89 / 17.98 | 27.43 / 27.97 | 109.7 / 111.9 s | 0.961 / 0.961 |
| A raw01 / raw02 | 16.35 / 14.64 | 23.43 / 24.20 | 93.7 / 96.8 s | 0.957 / 0.957 |
| B raw02 | 14.25 | 23.69 | 94.8 s | 0.981 |
| C raw01 / raw02 | 13.49 / 15.06 | 23.06 / 25.68 | 92.3 / 102.7 s | 0.976 / 0.981 |
| D raw01 / raw02 | 18.03 / 18.36 | 25.50 / 26.72 | 102.0 / 106.9 s | 0.983 / 0.976 |

- **Coverage:** at or above the old code's on every run.
- **Cost:** windows are wider, so more candidates per precursor in search; the effect on IDs is not measured.
- **Still open:** arm B raw01 still fails in stage 1, which has no RT information and is a different failure.
