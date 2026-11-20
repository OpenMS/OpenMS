# PLAN_IMNARROW: narrower calibrated 1/K0 window on the r ≈ 7 entrapment (frozen 20 Nov 2026)

**Question.** Does extracting within 0.6 × the calibrated 1/K0 half width raise IDs at equal true FDR on the large
library, without breaking own-q calibration? The aim is less background mixed into faint peaks, with decoys treated
by the same code.

**What was seen before this plan.** PRECHECK_IMNARROW, r ≈ 1.56 HeLa D, seed 0 only, `-im_tolerance_scale` sweep
(the calibrated half width there is ±0.052 1/K0 = 4 σ of the anchor residual, σ 0.0131):

| scale | N at FDP_c 1 % | lower / combined bound | pair win | wall |
|---|---|---|---|---|
| 1.0 (default, byte-identical) | 12,357 | 0.79 / 1.49 % | 0.524 | 5:03 |
| 0.8 | 12,806 | 0.69 / 1.31 % | 0.507 | 5:08 |
| 0.6 | 13,388 | 0.57 / 1.08 % | 0.501 | 6:13 |
| 0.5 | 13,666 | 0.60 / 1.13 % | 0.518 | 6:11 |
| 0.4 | 13,103 | 0.72 / 1.37 % | 0.510 | 7:38 |

- The faintest quartile gains most: at decoy FPR 1e-4, Q1 goes from 0.143 to 0.166 and Q2 from 0.262 to 0.298 (scale 0.6).
- **Choice of 0.6.** In run terms it is 2.4 σ of the calibration's own anchor residual (±0.031 here). It also matches
  the p95 of the 1/K0 error against DIA-NN in PRECHECK_IMPICK (0.031). That figure uses DIA-NN labels, so the choice
  is not label-blind. It was also made after seeing this sweep: 0.5 had more IDs on this seed, and 0.6 was taken as
  the value with the best calibration and pair win rather than the best count. Both are disclosed weaknesses of
  this freeze.
- No r ≈ 7 run with `-im_tolerance_scale` exists before this plan. The limits below were chosen knowing the r ≈ 1.56
  seed 0 result.

## Setup (fixed)
- **Binary:** OpenMS 0942ffbdf3, container `entlib-pdsif-20261118/opendia-pd10-0942ffbdf3.sif`. Its OpenDIAWorkflow
  binary has the same sha256 as the local build (5e09774d…). At scale 1 the output is byte-identical to the previous
  default on r ≈ 1.56 HeLa D, and that default equals PLAN_EXTRA NEW. 25/25 ctests pass.
- **Library:** the PLAN_ODP predicted library (sha256 76f2a6c2…), as `-lib`.
- **Arms:**
  - CUR is the current default.
  - NEW adds `-im_tolerance_scale 0.6`.
- **Runs:**
  - Seeds 0-4, D first (2737.d), 8 threads, 128 GB.
  - H (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py` (sha256 bb03387d…). Entrapment labels are read only there.

## Criteria (on D, pooled or median over the 5 seeds)
1. **Sensitivity:** NEW median N at FDP_c 1 % ≥ 1.03 × CUR median.
2. **Calibration not materially worse:**
   - NEW pooled lower bound at own q ≤ 0.01 ≤ CUR pooled lower bound + 0.0020.
   - NEW pooled combined bound ≤ CUR pooled combined bound + 0.0025.
3. **Pair fairness:** NEW P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.
4. **Cost:** NEW median wall ≤ 1.4 × CUR median wall, and peak RSS ≤ 64 GB.

## Outcome
- **All pass:** report H. Making 0.6 the default is then a separate commit, with ctests and a byte check.
- **Any fail:** report as is. The option stays at 1 and nothing is retuned.
