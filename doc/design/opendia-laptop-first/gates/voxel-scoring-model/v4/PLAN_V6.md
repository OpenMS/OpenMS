# Voxel v6: per-fold decoy-tail calibration of the final score (written before any v6 scores were computed)

**Context:**
- `RESULT_V4P2.md` limit 1: the 1 % cut is unstable on raw 02.
- `RESULT_SEQ.md`: DVS collapsed when the folds were pooled.
- The evaluator pools the 3 folds after a per-fold decoy median/MAD scaling, which matches the bulk of each fold's decoy winners but not the tail.
- Pre-plan look (this file's author, existing scores only): q99.9 of the scaled decoy winners differs between folds. Raw 02 D11FW gives 75 / 26 / 44, raw 01 37 / 39 / 37. Where one fold's tail is heavier, its decoys enter the pooled list early and move the 1 % cut.

## Transform (label-free; no training; same rule for targets and decoys)

For each fold k, let S_k(s) = (number of decoy precursor winners of fold k with score ≥ s + 1) / (number of decoy winners of fold k + 1), the empirical decoy survival.

- **Below t_k:** s' = −log S_k(s), linearly interpolated in s between decoy scores. t_k is the score at the top 1 % of fold k's decoy winners.
- **Above t_k:** s' = −log S_k(t_k) + (s − t_k) / β_k, with β_k the mean excess over t_k of those top-1 % decoys (exponential tail), so the transform is continuous and strictly increasing.
- **Scope:** applied to every window of the fold, so window choice within a precursor and ranking within a fold are unchanged. Windows set to −1e30 by the W rule stay at −1e30.

## Arms

`calib_tail.py`, no retraining:

| Arm | Source |
|---|---|
| D11FWT | D11FW + transform |
| D11FPWT | D11FPW + transform |

## Pass rule

Each arm, both raws; same as PLAN_V5.

1. PLAN_V2 rules 1–4 against BASE.
2. N@2 % above D11FW on both raws.
3. Null-winner AUC ≤ D11FW's + 0.005.

**Reported only:**
- N at 0.5 / 1 / 1.5 %;
- exact FDP_c;
- the funnel;
- the paired bootstrap against D11FW and D11.

An arm that passes becomes the base.
