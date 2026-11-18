# PLAN_PGENT2 design verdict: FAIL on criterion 1; the excess over 1 % is real, not only near-homologs

> **Correction (BOUNDS.md):** FDP_c here is the combined estimator, an *upper bound*. The lower bounds (0.87 % CUR, 0.97 % NEW) are below 1 %, so "q is about 1.5 times optimistic" is not shown. The true FDP lies between about 0.9 % and 1.6 %, which is inconclusive. The verdict on the frozen criterion stands.

Written 18 November 2026, after the D design runs and before any H run.
- **Library:** `hela_ath2.pqp` (sha256 e1ac12b8…), r = 1.560.
- **Evaluation:** `eval.py`, written to `EVAL_D.txt`.

| D, seeds 0-4 | N at q ≤ 0.01, median | entrapments, pooled | FDP_c at q ≤ 0.01, pooled | N@FDP_c 1 %, median (range) |
|---|---|---|---|---|
| CUR | 11,313 | 492 / 56,574 | 0.0143 | 10,576 (10,358-10,858) |
| NEW | 13,355 | 643 / 66,544 | **0.0159** | **12,388** (12,314-12,545), +17.1 % |

1. **Criterion 1 fails:** NEW FDP_c 0.0159 > 0.012.
2. **Criterion 2 passes:** +17.1 %. Every NEW seed is above every CUR seed.
3. Criteria 3 and 4 (Astral, restated) pass.

**Verdict: FAIL.** `-pg_scores` stays opt-in, and nothing is retuned. H is not run: the verdict cannot change, and the
plan only requires H to be reported for a pass or for its record. Running it would cost about 2 h on codon for a number
that does not decide anything. It can be run on request.

## Correction to PRECHECK_FOREIGN
- That precheck filtered near-homologs **after** the search, removing them from both N_e and r. It reported 1.16 % (CUR)
  and 0.97 % (NEW). This plan filters them **at construction**, and gets 1.43 % (CUR) and 1.59 % (NEW) over 5 seeds.
- So the post-hoc filter flattered both arms. Removing accepted near-homologs took out the entrapments most likely to
  score well. The remaining set was then compared with an r that still reflected weaker candidates.
- With near-homologs gone from the library itself, the remaining Arabidopsis entrapments still pass at about 1.4-1.6 %.
- The precheck's conclusion that "q is close to calibrated" is withdrawn. **OpenDIA's q ≤ 0.01 on diaPASEF with a
  predicted library gives about 1.4-1.6 % true FDR**, measured with a foreign, homolog-free entrapment. That is about
  1.5 times nominal.
- `-pg_scores` adds about 17 % IDs at equal true FDR, and is slightly less well calibrated than CUR: +0.16 percentage
  points, within one run's noise.

## What this means for option (b)
The calibration defect is real. It survives every benchmark variant tried: shuffled, foreign, and homolog-free. It is not
explained by:
- decoy type (pseudo-reverse versus shuffle);
- library properties (charge, m/z, length, RT, IM: AUC 0.49-0.51);
- parent interference.

The candidates still open are in the scoring and q pipeline itself:
1. The per-fold decoy-tail calibration (`calibrateTail`), whose exponential tail is fitted on the top 1 % of decoys.
2. Self-training positives at q 0.15 / 0.05, which on a low-signal predicted library include more false targets.
3. The "W rule" that keeps one window per precursor. A target picks its best of up to 5 windows, and so does a decoy.
   With multiple chances a random target and a random decoy should be symmetric, so this is unlikely, but it is not
   tested.

The next step is a run-derived diagnostic on the existing D outputs: rho per q band for each of those stages
(pre-calibration raw score, linear part, calibrated score). That locates where the asymmetry enters, without new runs.
