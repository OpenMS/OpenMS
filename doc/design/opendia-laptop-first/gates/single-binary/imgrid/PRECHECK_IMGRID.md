# PRECHECK_IMGRID (seed 0, HeLa design run D, 2026-11-17)

The change (uncommitted code, `OpenDIAWorkflow.imgrid`) is on top of 6dd893b954 (RTCENTRE). It has two parts.
- **12-fragment grid:** with a mobility calibration, the grid sums only points within the calibrated 1/K0 tolerance of
  each precursor (the same tolerance the search uses). This also changes the interference inputs.
- **Interference:** two precursors link only if their calibrated 1/K0 differ by at most 2 × the tolerance, that is,
  only if their mobility windows overlap.

Without a mobility calibration (Astral) nothing changes.

| D, seed 0 | N at q ≤ 0.01 | agreement with DIA-NN |
|---|---|---|
| RTCENTRE (6dd893b954) | 14,927 | 0.9335 |
| + IMGRID | 16,021 | 0.9307 |

Precursor-winner target-vs-decoy AUC by feature, RTCENTRE → IMGRID:

| feature | f12_support | f12_coelution | f12_libcos | i_libcos_resid | d_max |
|---|---|---|---|---|---|
| AUC | 0.593 → 0.624 | 0.621 → 0.661 | 0.685 → 0.714 | 0.615 → 0.652 | 0.417 → 0.401 |

Every grid and interference feature separates targets from decoys better once the extraction respects mobility.
Agreement with DIA-NN drops by 0.003, which is why PLAN_IMGRID keeps an agreement rule.
