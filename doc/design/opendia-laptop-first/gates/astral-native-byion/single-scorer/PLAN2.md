# Single TOPAZ scorer, round 2: augmentation against overfitting

Written 2026-11-09 17:45 UTC, after `RESULT.md` (raw 01 S0 FAIL) and before any round 2 run. The code and binary are unchanged (`33d687b330`, `binary.sha256`). Raw 03 stays untouched.

## Why

In S0 round 0, the about 600 seeds per fold were trivially separable with the heuristics (validation loss about 2e-5, seed AUC 1.0), and only 347 targets passed q ≤ 0.01. The score iterations made this worse (149), then crashed.

Augmentation should make the training task harder and closer to the real one:
- weak targets (s × target + decoy window, s in [0.1, 1]);
- missing non-trace features (dropout).

## Arms (sequential, about 16 min and 25 GB each)

Common: the `run.sh` command with `-ss_iterations 0`. There are no score iterations, so q comes from the round-0 TOPAZ logit.

| Arm | Raw | Variable |
|---|---|---|
| **A1** | raw 01, then raw 02 | `-augment_probability 0.5 -augment_scale_min 0.1` |
| **C0** | raw 01 only, run last | no augmentation. This is the control for A1 and gives the N@1% that S0 never wrote. |

## Decision rule

The rule is the same as `PLAN.md` rules 1–4, unchanged and applied to A1 on both raws:
1. N@1% ≥ 1050 (raw 01) and ≥ 1324 (raw 02).
2. FDP_c at q ≤ 0.01 is ≤ 0.015.
3. FDP_c at q ≤ 0.05 is ≤ 0.070 (raw 01) and ≤ 0.079 (raw 02).
4. The run completes, with peak RSS ≤ 30 GB and wall time ≤ 35 min.

Outcomes:
- If A1 passes, it is the candidate replacement for CNN + LDA. It is not the default without the raw 03 gate, which needs approval.
- **Augmentation effect (reported, not a pass condition):** A1 against C0 on raw 01, for N@1% and q ≤ 0.01 counts. Augmentation counts as helping if A1's N@1% is at least 1.05 × C0's N@1%.
- If A1 fails but helps, the next step is a stronger augmentation (a higher p or a lower s_min), declared separately.

Reported, as in `PLAN.md`: N at FDP_c 0.5 / 1.5 / 2%, the q counts, R per fold, validation loss, `augmented_first_epoch`, time and RSS. All q-derived counts are decoy-count estimates.
