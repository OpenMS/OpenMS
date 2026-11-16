# PLAN_RTPROFILE: recalibrate from the run's own IDs, with an RT-dependent window

Written 16 November 2026, before the implementation and before any run. The evidence is in `PRECHECK_RTPROFILE.md`.
The user asked for (a) a window whose width follows the gradient and (b) more calibration points at the gradient ends,
taken from the target IDs found there.

## Why the calibration has no anchors at the ends (raw01 calibration diagnostics)
The panel already covers the ends: it samples targets evenly over library RT, 17,182 at library RT < 0.03 and 55,244 at
0.85-0.88.
- Calibration stage 3 still accepts almost none of them at the anchor FDR: 0 / 1 / 2 anchors in the end regions, and
  only 14 / 15 / 34 targets even at q ≤ 0.2.
- Outside the first-pass map's support, stage 3 searches to the run edge with its co-elution score alone, which cannot
  separate them from decoys.
- The main analysis, with its 24 features and network, does identify them: raw01 has 71 IDs at q ≤ 0.01 below library
  RT 0.07, and 89 above 0.8.
- So the end anchors should come from the analysis IDs, not from the calibration score.

Offline check (raw01 seed 0): refitting the map on the pass-1 IDs, with 10-fold CV, gave these 95th-percentile errors
before the first anchor.

| run | before (s) | after (s) |
|---|---|---|
| raw01 | 130 | 72 |
| raw03 | 36 | 22 |
| raw02 | 122 | 76 (median worse: 36 → 48) |

The interior is unchanged.

## Design (option `-id_recalibration` in OpenDIAWorkflow)
1. **Pass 1.** Run as now: calibration, then search with window ±T, population, features, scorer.
2. **ID anchors.** Take the pass-1 target precursors at run q ≤ 0.01: their apex and their library RT.
   - Pass 2 runs only if every fold has at least 20 ID anchors outside it. Otherwise pass 1 is reported.
3. **Cross-fitted maps.** For each fold k, M_k = `fitMonotoneRTMap` on the ID anchors of the **other two folds**.
   - A precursor (and its decoy, which shares the fold) is centred with its fold's map.
   - So no target's own pass-1 peak sets its own centre.
4. **Width profile, cross-fitted.**
   - Held-out residuals: e_i = |apex_i − M_{fold(i)}(library_rt_i)| for every ID anchor, against the map that did not
     see it.
   - For fold k, using only the anchors of the other folds:
     - σ_k(x) = 1.4826 × median of e over the round(√n) anchors nearest in library RT, where n is that anchor count.
     - σ̄_k = 1.4826 × median of all of them.
   - Half width w_k(x) = T × σ_k(x) / σ̄_k.
   - T stays the calibration's selection-free level, so the step changes the **centre and the shape** of the window,
     not its overall level. It narrows where the error is small and widens where it is large: the ends and sparse
     regions.
5. **Pass 2.** Search again with per-precursor centres and widths, then population, extraction, features and scorer as in
   pass 1.
   - The scorer seeds stay the calibration anchors, as before.
   - The RT error feature is computed against the new centre.
6. **Output.**
   - Pass 2 → `<run>.precursors.tsv`; pass 1 → `<run>.precursors.pass1.tsv`.
   - The log reports the anchor counts and the w_k range.

Only target/decoy status and q-values are used, never entrapment labels. The only constants are the existing ones:
the q ≤ 0.01 report level, the 20-anchor minimum of the map fit, √n as in the map's end-slope rule, and T.

## Arms, runs, criteria (as PLAN_RTWIN2)
- **Arms:** CUR = pass 1, NEW = pass 2, from the same run. Self-calibrating, arm-E library `library.pqp` and
  `library12_full.pqp`, seeds 0-4. Codon, 8 threads.
- **Runs:** raws 01 and 02 are the design runs. Raw 03 is held out, evaluated once after the verdict is committed.
- **Evaluation:** `eval_run.py`, with entrapment labels used only there.
- **Criteria,** on raws 01 and 02 separately:
  1. Median N@FDP_c 1 %: NEW ≥ CUR − 1 %.
  2. FDP_c at q ≤ 0.01 ≤ 0.010 for every NEW seed.
  3. Additionally, NEW ≥ CUR + 2 % on at least one raw.
- **Outcome:** pass makes `-id_recalibration` the default. Fail removes the code. Nothing is retuned.
