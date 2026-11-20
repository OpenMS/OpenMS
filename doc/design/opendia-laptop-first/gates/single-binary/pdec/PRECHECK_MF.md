# PRECHECK_MF: fragment-level, mobility-aware matched filter for faint right peaks (20 Nov 2026, diagnostic)

**Question.** The missing DIA-NN IDs are mostly right peaks that score too low (B1, 62.5 % of the missing). Can per-fragment matched filtering, meaning robust per-fragment S/N, leave-one-out shape agreement and fragment dominance on mobility- and ppm-weighted traces, separate them from decoys better than the 36 default features?

## Setup

- **Sample:** faint2, which has 6,000 DIA-NN non-entrapment q ≤ 0.01 targets and 30,000 decoys on r ≈ 1.56 HeLa D. Raw points come from `maps.xipm` (±15 ppm, ±0.08 1/K0).
- **Baseline:** the current default 36 features (10 ppm trace, `-im_tolerance_scale 0.6`), from `faintcpp-20261120/runs/D_cur_s0`.
- **Script:** `precheck_mf.py` (sha256 `48ba322d…`), codon job 39077071. The gate was written into the script header before any result. Amendment 1 added the tensor-feasible `mt_*` variant after job 39074453 was cancelled at 1:15, before any feature or result existed.
- **Variants:**
  - `mf_*` uses raw points on the native frame grid (about 0.95 s), weights each point by its calibrated ppm error and by a Gaussian (σ 0.008) around the window's own consensus mobility apex.
  - `mt_*` computes the same statistics from what the workflow's trace tensor already holds: 10 ppm, 32 RT bins of 2.2 s, 8 IM bins, and IM weighting on bin centres around the tensor's IM apex.
- **Arms:** a GBM cross-fitted on the tool's folds, 2 seeds.
  - **Ceiling:** DIA-NN right windows against decoys.
  - **Tool proxy:** label-free, with positives being the tool's own q ≤ 0.01 targets.
- **Metric:** TPR of DIA-NN right windows at decoy FPR 3e-3, over all IDs and over the 2,537 DIA-NN IDs the tool misses.
- **Gate:** +0.030 in the ceiling arm and +0.020 on the missed subset in the tool arm, on both seeds. core5 had +0.022 and recovered no IDs, so the bar sits above it.

## Result

| Arm | Ceiling, all (s0 / s1) | Ceiling, missed | Tool, all | Tool, missed |
|---|---|---|---|---|
| 36 | 0.736 / 0.731 | 0.268 / 0.255 | 0.719 / 0.716 | 0.225 / 0.217 |
| +core5 | +0.013 / +0.016 | +0.038 / +0.046 | +0.013 / +0.006 | +0.025 / +0.013 |
| **+mf** | **+0.066 / +0.071** | **+0.195 / +0.208** | **+0.053 / +0.038** | **+0.143 / +0.113** |
| **+mt** | **+0.034 / +0.037** | **+0.102 / +0.111** | **+0.038 / +0.033** | **+0.103 / +0.091** |
| +mf+core5 | +0.070 / +0.072 | +0.207 / +0.212 | +0.049 / +0.041 | +0.137 / +0.123 |

**Both variants pass the gate.** This is the first feature block that moves the label-free tool proxy by more than about 0.01, and it does so mainly on the missed IDs. On the missed subset, the tool-proxy TPR at decoy FPR 3e-3 roughly doubles with `mf`, from 0.22 to 0.33-0.37.

Single-feature AUCs on missed right windows against decoys:

| Feature | mf | mt |
|---|---|---|
| llr | 0.830 | 0.819 |
| n_z2 | 0.812 | 0.816 |
| z_med | 0.812 | 0.817 |
| loo_med | 0.760 | 0.766 |
| cos | 0.736 | 0.723 |
| dom | 0.300 | 0.344 |

For comparison, the best existing features on B1 were about 0.70-0.73.

## Reading

- **The signal is fragment-level robustness, not mobility alone.** core5 (mobility only) adds little. The gain comes from asking how many fragments individually rise above their own local background (`n_z2`, `z_med`, `llr`), whether each fragment's shape agrees with the others (`loo`), and whether one fragment dominates (`dom`). Pooled sums let one interfering fragment carry a window. Per-fragment robust statistics do not.
- **mt carries about half of mf's ceiling gain and most of its tool-proxy gain.** The remainder needs the native frame grid (2.3× finer RT) and ppm weighting, which require a new raw-point grid in the extraction pass.
- **Caveat:** these are GBM ceilings and a GBM label-free proxy on one run. The workflow scorer is an MLP, and earlier ceiling gains realised less. The decision rests on the r ≈ 1.56 workflow check, then on a frozen r ≈ 7 plan.

## Decision

As amended before results: implement `mt` first (`-matched_filter`, opt-in, no extraction change), then run a two-seed r ≈ 1.56 workflow check against the current default. `mf` is the follow-up if `mt` realises a gain, because it needs the extraction pass to keep a finer, ppm-weighted grid.
