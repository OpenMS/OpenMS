# Single TOPAZ scorer (no CNN + stage-2 LDA): raws 01 and 02

Written 2026-11-09 17:05 UTC, before any single-scorer run. Code: `33d687b330` (`-topaz_features`, `-augment_probability`).
Raw 03 stays untouched. DIA-NN and OSW are not used for training. They are only reference counts.

## Question

Can one learner replace grid2d CNN + stage-2 LDA without losing IDs or calibration? The learner is TOPAZ width 8 with the heuristic groups in its head, semi-supervised score iterations, and q from its own logit.

## Arms (identical except for the variable named)

Common to both arms, matching the frozen native byion command except where noted:
- `-cnn_arch topaz -topaz_width 8 -weight_decay 0.01`. This is the compact layout of the earlier arm C.
- `-topaz_features search+boundary+byion` replaces `-stage2 -stage2_feature_set search+boundary+byion`.
- `-ss_iterations 2 -ss_fdr 0.01`: PyProphet-style retraining on target winners at q ≤ 0.01 against decoy winners.
- `-paired_tdc -pooled_tdc`, with the same calibration, seeds, windows, search table, folds 3, seed 42, epochs 40, patience 6 and batch 32.

| Arm | Variable |
|---|---|
| **S0** | no augmentation |
| **S1** | `-augment_probability 0.5 -augment_scale_min 0.1` |

Order: raw 01 S0, raw 01 S1, raw 02 S0, raw 02 S1, run one at a time because each peaks at about 25 GB.

Before the arms run, the HeLa default command (`cnn-topaz-design/hela-default-cmp4/run.sh`) must give a byte-identical `crossfit.tsv` with the new binary. If it does not, the chain stops.

## Metric

The same as the native byion reading (`astral-raw02-byion-design-20261109/evaluate.py`).
- **N@1%:** the number of in-window target precursor winners at entrapment FDP_c ≤ 1%, ranked by the score. Here FDP_c = N_e(1 + 1/r)/N, with r = 148983/151239.
- **Pooled paired q:** the tool's TDC on the final-round logit. From it, FDP_c at q ≤ 0.01 and at q ≤ 0.05.

Base (the native byion design, CNN + LDA):

| | raw 01 | raw 02 |
|---|---|---|
| N@1% | 1029 | 1298 |
| FDP_c at q ≤ 0.01 | 0.013 | 0.007 |
| FDP_c at q ≤ 0.05 | 0.065 | 0.074 |

OSW + Percolator reference N@1%: 1151 (raw 01) and 1287 (raw 02).

## Decision rule (fixed now)

An arm **passes** only if all of these hold on **both** raws:
1. N@1% ≥ 1.02 × base, that is ≥ 1050 on raw 01 and ≥ 1324 on raw 02.
2. FDP_c at q ≤ 0.01 is ≤ 0.015.
3. FDP_c at q ≤ 0.05 is ≤ base + 0.005, that is ≤ 0.070 on raw 01 and ≤ 0.079 on raw 02. This protects raw 02 calibration.
4. The run completes, with peak RSS ≤ 30 GB and wall time ≤ 35 min.

Outcomes:
- If both arms pass, the one with the higher summed N@1% is preferred. On a tie within 1%, S0 is preferred because it is simpler.
- A passing arm becomes the candidate replacement for CNN + LDA. It does **not** become the default: that needs the raw 03 gate, which needs the user's approval.
- If neither arm passes, CNN + LDA stays the native design and the new options stay opt-in and experimental.

Reported but not part of the decision:
- N at FDP_c 0.5 / 1.5 / 2%;
- q ≤ 0.01 / 0.05 counts;
- per-round `target_precursors_q01`;
- the null-symmetry ratio R at the decoy-winner median per fold (as in `cnn-topaz-design/PLAN.md`);
- per-fold validation loss and `augmented_first_epoch`;
- wall time and RSS.

All q-derived counts are decoy-count estimates.
