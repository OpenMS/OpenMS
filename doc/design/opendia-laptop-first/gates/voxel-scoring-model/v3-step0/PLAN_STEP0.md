# Voxel v3, step 0: two checks before any new training (written before the checks run)

No training of the voxel model. The verdict of PLAN_V2 does not change. Inputs: the frozen v2 scores (`scores_V0.npy`, `scores_FS.npy`), the v2 voxels, and the base features (`features.tsv`, with `cnn_z` and the 11 heuristic features). Raws 01 and 02.

## Check A: does the voxel score depend on brightness?

- Brightness of a window = log(1 + sum of the 6 fragment traces over the 32 RT bins), from `xic.npy`.
- Positives: the calibration anchor windows (group `anchor` in `train_windows.tsv`). Every score is out-of-fold for them.
- Negatives, two sets, both from all decoy windows of the population (not only winners):
  - R: 20 random decoy windows per anchor;
  - M: 20 decoy windows per anchor, drawn from the same brightness bin (20 quantile bins of the anchor brightness; decoys outside the anchor range are not used).
- Scores: BASE (`stage2_score`), `cnn_z`, V0, FS. For each, AUC(anchor > R) and AUC(anchor > M). Drop = AUC_R − AUC_M.
- Check: brightness alone must give AUC near 0.5 on M (shows that the matching works).
- **Shortcut confirmed** for an arm if, on both raws, its drop is ≥ 0.01 and ≥ 2 × the BASE drop.
- A second, weaker view (it is biased to the base): target winners at base q ≤ 0.01 that are not anchors, against M-matched decoy windows.

## Check B: does the voxel score add information to the base?

- Python copy of the base stage-2 LDA (`crossfitStage2LDA`): base folds (`fold` column), standardise on the training rows, 3 rounds, round 1 ranks by feature 0 (`cnn_z`) at q ≤ 0.15, then q ≤ 0.05, positives = target winners, negatives = all decoy winners of the training rows, ridge 1e-3. Population = the v2 population (pairs with max transition support ≥ 3).
- Arms:
  - B0: the 12 base features (`cnn_z` + 11). This is the comparison base of check B.
  - B1: B0 + V0 score.
  - B2: B0 + FS score.
- Note: the V0 / FS scores are out-of-fold under the v2 folds, `cnn_z` under the base folds. Both are honest out-of-fold scores, so a stacked LDA can use them.
- Metrics: as `evaluate_v2.py` (pooled paired TDC, N at FDP_c 0.5/1/1.5/2 %, q ≤ 0.01 / 0.05 N and FDP_c, entrapment AUC), with a paired bootstrap of B1 − B0 and B2 − B0 (1000 resamples, seed 5).
- Reproduction check: B0 N@1 % within ± 5 % of the population BASE (1024 / 1288). If not, report it and read B1 / B2 only against B0.
- **The voxel score adds information** if, for B1 on both raws: the 90 % interval of ΔN@2 % is above 0, the ΔN@1 % interval is not entirely below 0, FDP_c at q ≤ 0.01 ≤ 0.015, and the entrapment AUC ≤ the B0 AUC + 0.005.

## Use of the result

- A confirmed (and B fails): v3 steps 1 and 2 go ahead (remove the brightness shortcut, matched and symmetric negatives).
- B passes: a separate result. It is a candidate C++ change (voxel score as one LDA feature). It needs its own plan and the user's approval.
- A not confirmed and B fails: report to the user before v3. The brightness cause is then not supported.
