# PLAN_IM2D: RT × ion-mobility candidate picking on the r ≈ 7 entrapment (frozen 19 Nov 2026)

**Question.** Does `-im_peaks 8` (OpenMS 07131138b8) raise IDs at equal true FDR on the large library, where co-eluting
precursors at other mobilities are 4.5× denser than at r ≈ 1.56, without breaking own-q calibration?

**What was seen before this plan.** PRECHECK_IM2D, r ≈ 1.56 HeLa D, seed 0 only:
- right-RT candidates for DIA-NN IDs went from 84.1 % to 93.5 %;
- N at FDP_c 1 % went from 12,357 to 12,788 (+3.5 %);
- the lower / combined bounds went from 0.79 / 1.49 % to 0.70 / 1.32 %;
- wall time rose 23 %.

No r ≈ 7 run with `-im_peaks` exists before this plan. The limits below were chosen knowing the r ≈ 1.56 seed 0 result.

## Setup (fixed)
- **Binary:** OpenMS 07131138b8, container `entlib-pdsif-20261118/opendia-pd9-07131138b8.sif` (sha256 312a8eb3…).
  - The default (CUR) arm of this binary is byte-identical to 969b474a51 on r ≈ 1.56 HeLa D. That is the evidence-block
    default, which equals PLAN_EXTRA NEW byte for byte (checked on r ≈ 1.56).
- **Library:** the PLAN_ODP predicted library (sha256 76f2a6c2…), as `-lib`.
- **Arms:**
  - CUR is the current default.
  - NEW adds `-im_peaks 8`, keeping the other im options at their defaults (`-im_extra_peaks 2`, `-im_min_support 1`).
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
4. **Cost:** NEW median wall ≤ 1.5 × CUR median wall, and peak RSS ≤ 64 GB.

## Outcome
- **All pass:** report H. Making `-im_peaks 8` the default is then a separate commit, with ctests and a byte check.
- **Any fail:** report as is. The option stays opt-in and nothing is retuned.
