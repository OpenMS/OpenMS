# Voxel v3 phase 1: is the training procedure the bottleneck? (written before any run)

Context: one single model is the target (no second-stage LDA). In v2, F0 (voxel + 11 heuristics, one network) gave 697 / 698 at FDP_c 1 %. A ridge LDA on the same windows with 12 base features gave 1181 / 1345 (step 0, B0). A network that gets the same inputs as a linear model must be able to match it. If it does not, the cause is the training procedure, and a better voxel input cannot fix it.

Population, folds, T0 rows, metrics and bootstrap: as PLAN_V2 (`evaluate_v2.py`, against the population BASE). Inputs: only the 11 heuristic features (no `cnn_z`, no voxel), standardised on the training folds. Entrapment labels are not used in training.

## Arms (one variable at a time)

| Arm | Model | Schedule | Negatives in rounds 1–2 |
|---|---|---|---|
| L11 | ridge LDA (1e-3) | round 0 on T0, rounds 1–2 at q ≤ 0.15 / 0.05 (as v2) | all decoy winners of the training folds |
| N11 | the v2 head network (11 → 128 → 64 → 1, dropout 0.1) | as v2 (round 0 on T0, 2 rounds, early stopping on family-split validation) | as v2: max(6000, 4 × positives) sampled decoy winners |
| N11a | as N11 | as N11 | all decoy winners (class-balanced weights, as `fit()` already does) |

## Reading (no pass to C++ from this phase)

- N11 vs L11: a network against a linear model with the same inputs and schedule.
- N11a vs N11: effect of the negative subsampling.
- The **procedure is a bottleneck** if L11 − N11 at N@1 % is > 5 % of L11 on both raws. Then phase 2 uses the better procedure (N11a if it closes the gap to ≤ 5 %).
- If N11 is within 5 % of L11, the procedure is not the cause. Phase 2 then works on the voxel representation only.
- Phase 2 (voxel encoder with label-free pretraining on all windows, then the chosen procedure) needs its own frozen plan.
