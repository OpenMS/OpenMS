# PLAN_EXTRA: `-extra_scores` on the r ≈ 7 entrapment (frozen 19 Nov 2026)

**Question.** Does the experimental evidence block (`-extra_scores`, 12 features: core-window fragment and MS1 mass
error, MS1 co-elution, isotope and area, MS2 area, ion-mobility centroid, co-elution and central share) raise IDs at
equal true FDR on a large library without breaking own-q calibration?

**What was seen before this plan.**
- PRECHECK_IDGAP (r ≈ 7 D seed 0, existing features only): a GBM ceiling trained on DIA-NN labels gains nothing over
  the tool score, so the gap is in the features.
- PRECHECK_EVIDENCE, r ≈ 1.56 HeLa D seed 0 only (OpenDIA-predicted decoys):

  | arm | N at q ≤ 0.01 | N at FDP_c 1 % | lower / combined bound |
  |---|---|---|---|
  | CUR | 9,953 | 9,863 | 0.55 % / 1.04 % |
  | NEW | 13,034 | 12,357 | 0.79 % / 1.49 % |

  Pair win is 0.505 for CUR and 0.524 for NEW. No r ≈ 7 run with `-extra_scores` exists before this plan. The limits
  below were chosen knowing the r ≈ 1.56 seed 0 result.

## Setup (fixed)
- **Binary:** OpenMS 83a0f6c14a, container `entlib-pdsif-20261118/opendia-pd8-83a0f6c14a.sif` (sha256 c6bcac1b…).
  - For the default arm the binary is byte-identical in output to PLAN_ODP CUR. This was checked through pd7 on r ≈ 7
    D/H, and on plain HeLa D for 83a0f6c14a.
- **Library:** the predicted library that PLAN_ODP wrote, `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`
  (sha256 76f2a6c2…). It is passed as `-lib` without `-decoys predict`. Reuse is output-identical (DIAG_SCALE, `-decoy_cache`).
- **Arms:** CUR is the default and NEW adds `-extra_scores`.
- **Runs:** seeds 0-4, D first (20200505 … 2737.d), 8 threads.
- **H** (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py` (sha256 bb03387d…) on the tool's own output. Entrapment labels are read only there.

## Criteria (on D, pooled or median over the 5 seeds)
1. **Sensitivity:** NEW median N at FDP_c 1 % ≥ 1.10 × CUR median.
2. **Calibration not materially worse:** NEW pooled lower bound at own q ≤ 0.01 ≤ CUR pooled lower bound + 0.0025, and
   NEW pooled combined bound ≤ CUR pooled combined bound + 0.0030.
3. **Pair fairness:** NEW P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.

## Outcome
- **All pass:** report H, then propose `-extra_scores` as the default in a separate commit with ctests. Absolute
  calibration (PLAN_ODP limits) stays an open item for both arms.
- **Any fail:** report as is. The flag stays opt-in and nothing is retuned.
