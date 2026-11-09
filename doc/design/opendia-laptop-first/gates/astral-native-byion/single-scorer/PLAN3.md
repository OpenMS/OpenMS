# Single TOPAZ scorer, round 3: decoy negatives

Written 2026-11-09 18:44 UTC, after `RESULT.md` rounds 1–2 and before any round 3 run. Code: `b75ad73a9c` (`-topaz_decoy_negatives`). Raw 03 stays untouched.

## Why

Every earlier arm (S0, A1, C0) fit its about 600 seeds per fold perfectly (held-out seed AUC 1.0). The held-out decoy tail then set q: C0 had 13 at N@1%, A1 had 0 at q ≤ 0.01.

The LDA, which passes, trains against about 177k decoy precursor winners per fold. TOPAZ trains against about 300 decoy seeds. The hypothesis is that the deep scorer needs the decoy distribution itself as negatives. Unlabeled decoy windows carry no target evidence, so this adds no label information about targets.

## Arm

**N1** = C0 (the PLAN2 command: TOPAZ width 8, weight decay 0.01, `-topaz_features search+boundary+byion`, `-ss_iterations 0`, no augmentation) + `-topaz_decoy_negatives 20000`.
- Each epoch draws 20000 random unlabeled training-family decoy windows as negatives.
- Validation adds 20000 fixed validation-family decoys.
- BCE pos_weight balances the classes.
- Order: HeLa default equivalence (byte-identical `crossfit.tsv` required, or the chain stops), then raw 01 N1, then raw 02 N1.

## Decision rule

The rule is the same as `PLAN.md` rules 1–4, unchanged, applied on both raws:
1. N@1% ≥ 1050 (raw 01) and ≥ 1324 (raw 02).
2. FDP_c at q ≤ 0.01 is ≤ 0.015.
3. FDP_c at q ≤ 0.05 is ≤ 0.070 (raw 01) and ≤ 0.079 (raw 02).
4. The run completes, with peak RSS ≤ 30 GB and wall time ≤ 35 min.

Outcomes:
- If N1 passes, it is the candidate replacement for CNN + LDA. It is not the default without the raw 03 gate, which needs approval.
- **Effect (reported, not a pass condition):** N1 against C0 on raw 01, for N@1% and q ≤ 0.01 counts.
- If N1 fails but improves q ≤ 0.01 by ≥ 2× over C0 (≥ 694), the hypothesis is supported. Next would be adding augmentation, or a larger n, declared separately.

Reported, as before: N at FDP_c 0.5 / 1.5 / 2%, R per fold, validation loss, time and RSS. All q-derived counts are decoy-count estimates.
