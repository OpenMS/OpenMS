# PLAN_CALSPREAD: selection-free RT spread with a member-shifted background

Written 15 November 2026, before any run of the design below.

## Problem
`OpenDIACalibration` stage 4 fits the RT tolerance (4 × `rt.spread_s`) as a mixture over the stage-3 target winners:
a Gaussian for true peaks plus a background. The two parts use different coordinates:
- the target residuals are taken to the final map;
- the decoy background is taken to the stage-3 window centre, which is the stage-1 or stage-2 map.

When the two maps differ (median 11.4 s on arm E raw01), the background edges do not line up with the targets' and the
Gaussian absorbs the mismatch. Measured against the robust spread of DIA-NN 1 % IDs around the same map:
- arm E raw01: 37.7 s against 27.2 s;
- arm E raw02: 30.9 s against 25.4 s (`single-binary/workflow/RESULT_WORKFLOW.md`).

Two offline variants were tried before this plan (`SPREAD_FIX*.txt`, descriptive):
- decoys taken to the final map;
- window frame without excluding anchors, which collapsed to 0-5 s on 2 of 10 runs.

They are not this design.

## Design (one variant, no tuning)
- **Data:** each stage-3 target winner m with its residual x_m to its own window centre, x_m in [-hw, hw], and the shift
  s_m = final map − window centre at the target's library RT.
- **Background:** the decoy winners' residual histogram b(x) in window coordinates, the same density used now. In
  window coordinates every member's chance peak has this density, so its edges match.
- **True peaks:** they sit at s_m + N(mu, sigma) in window coordinates.
- **Mixture:** p(x_m) = pi N(x_m − s_m; mu, sigma) + (1 − pi) b(x_m), fitted by EM.
- **Anchors excluded:** the stage-3 anchors, the targets accepted at the anchor FDR, are left out of the fit. They chose
  the final map, so their residuals to it are selected, and they are the spike that collapsed the earlier window-frame
  variant.
- **Start and stopping:** EM starts at sigma = stage-1 sigma, mu = 0, pi = 0.3, and runs to convergence.
- **Acceptance and fallback:** the fit is accepted under the current rules (≥ 20 peaks' worth, sigma < hw); otherwise
  the map sigma is used, as now.
- **Tolerance:** 4 × sigma, unchanged.

## Evaluation (labels and DIA-NN only here)
Offline first, on the 10 existing calibration outputs with `diag.tsv`: v2 arms A, C, D (raw01 and raw02), B raw02,
pairlib raw01 and raw02, and arm E raw01. Then in C++.

1. **Spread:** sigma within ±15 % of DIA-NN's robust spread around the same map on every run.
2. **Coverage:** the share of DIA-NN 1 % IDs within 4 sigma of the map is ≥ 0.95 on every run.
3. **End to end:** C++ and the self-calibrating OpenDIAWorkflow on both raws (arm-E library, seed 0).
   - FDP_c at q ≤ 0.01 must be ≤ 0.010 on both raws.
   - N@FDP_c 1 % must be at least the self-calibrated value of the current code (1718 / 1941).

The design passes if 1-3 hold. If 1 or 2 fail offline, the C++ change is not made, and the result is recorded as it is.
