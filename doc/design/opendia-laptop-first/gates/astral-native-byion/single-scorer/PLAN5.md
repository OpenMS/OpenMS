# Single TOPAZ scorer, round 5: target-window negatives

Written 2026-11-09 20:37 UTC, after `RESULT.md` rounds 1–4 and the diagnostics, and before any round 5 run. Code: `641a315bd6` (`-topaz_negative_source`). Raw 03 stays untouched.

## Why

Every TOPAZ arm ranks entrapment (false) target winners above decoy winners (AUC 0.52 per fold). The stage-2 LDA is at 0.50.

In N1 the extra negatives were all decoys, so "is a decoy" was a valid shortcut to the negative class. If the extra negatives are instead unlabeled target windows (mostly absent peptides and non-winning candidates), decoy-specific features no longer separate the classes. The scorer then has to learn peak evidence.

There is a cost. About 0.1% of target windows are true peaks labelled negative (about 6k true IDs out of about 5.4M windows), which is label noise.

## Arms (sequential, about 20 min and 25 GB each)

Both arms are N1 (the PLAN3 command: TOPAZ width 8, weight decay 0.01, `-topaz_features search+boundary+byion`, `-ss_iterations 0`, `-topaz_decoy_negatives 20000`, MLP head) plus:

| Arm | Variable |
|---|---|
| **T1** | `-topaz_negative_source target` |
| **M1** | `-topaz_negative_source mixed` |

Order: HeLa default equivalence (byte-identical, or the chain stops), then raw 01 T1, raw 01 M1, raw 02 T1, raw 02 M1.

## Decision rule

The rule is the same as `PLAN.md` rules 1–4, unchanged, applied per arm on both raws:
1. N@1% ≥ 1050 (raw 01) and ≥ 1324 (raw 02).
2. FDP_c at q ≤ 0.01 is ≤ 0.015.
3. FDP_c at q ≤ 0.05 is ≤ 0.070 (raw 01) and ≤ 0.079 (raw 02).
4. The run completes, with peak RSS ≤ 30 GB and wall time ≤ 35 min.

Outcomes:
- If both arms pass, the one with the higher summed N@1% is preferred.
- A passing arm is the candidate replacement for CNN + LDA. It is not the default without the raw 03 gate, which needs approval.
- **Mechanism check (reported, not a pass condition):** the entrapment-vs-decoy winner AUC (`decoyness.py`) per fold. The hypothesis is supported if T1's mean AUC is ≤ 0.51 on both raws (N1: 0.519 / 0.521).

Reported, as before: N at FDP_c 0.5 / 1.5 / 2%, q counts, R per fold, held-out seed AUC, time and RSS. All q-derived counts are decoy-count estimates.
