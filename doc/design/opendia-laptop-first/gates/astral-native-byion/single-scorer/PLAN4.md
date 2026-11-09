# Single TOPAZ scorer, round 4: linear head

Written 2026-11-09 19:32 UTC, after `RESULT.md` round 3 and before any round 4 run. Code: `f0f7a6205a` (`-topaz_linear_head`). Raw 03 stays untouched.

## Why

N1 (decoy negatives) ranked close to the LDA base: N@1% was 950 / 893 against 1029 / 1298. q from its logit was miscalibrated both ways: raw 01 q ≤ 0.01 had FDP_c 0.000 with 340 IDs, raw 01 q ≤ 0.05 had FDP_c 0.091, and seed AUC was still 1.0.

The calibrated reference (stage-2 LDA) is linear over 12 features. The hypothesis is that a linear head over the TOPAZ embeddings, co-elution statistics and heuristics overfits the seeds less and gives decoy-like scores to absent targets.

## Arm

**L1** = N1 (the PLAN3 command with `-topaz_decoy_negatives 20000`) + `-topaz_linear_head`.

Order: HeLa default equivalence (byte-identical, or the chain stops), then raw 01 L1, then raw 02 L1.

## Decision rule

The rule is the same as `PLAN.md` rules 1–4, unchanged, applied on both raws:
1. N@1% ≥ 1050 (raw 01) and ≥ 1324 (raw 02).
2. FDP_c at q ≤ 0.01 is ≤ 0.015.
3. FDP_c at q ≤ 0.05 is ≤ 0.070 (raw 01) and ≤ 0.079 (raw 02).
4. The run completes, with peak RSS ≤ 30 GB and wall time ≤ 35 min.

Outcomes:
- If L1 passes, it is the candidate replacement for CNN + LDA. It is not the default without the raw 03 gate, which needs approval.
- **Effect (reported, not a pass condition):** L1 against N1 on both raws. Calibration counts as improved if |FDP_c at q ≤ 0.05 − 0.05| shrinks on both raws.

Reported, as before: N at FDP_c 0.5 / 1.5 / 2%, q counts, R per fold, held-out seed AUC, time and RSS. All q-derived counts are decoy-count estimates.
