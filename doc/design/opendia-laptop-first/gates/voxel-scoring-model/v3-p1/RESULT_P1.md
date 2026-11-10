# Voxel v3 phase 1 result: the training procedure is a bottleneck

Plan: `PLAN_V3_P1.md` (sha256 `31f288fc…`, written 2026-11-10T14:39:03Z, before the run). The sha was checked after the run. Output: `EVALUATION_P1.txt`. Inputs: the 11 heuristic features only (no CNN, no voxel). Population, folds, schedule and metrics as v2.

| Raw | Arm | N@1 % | N@2 % | Δ@1 % vs BASE | Δ@2 % vs BASE | FDP_c q ≤ 0.01 | Entrapment AUC |
|---|---|---|---|---|---|---|---|
| 01 | BASE (CNN + LDA) | 1024 | 1314 | | | 0.013 | 0.5313 |
| 01 | L11 (LDA, 11 features) | **1267** | **1484** | +31..+368 | +68..+266 | 0.010 | 0.5310 |
| 01 | N11 (network, sampled negatives) | 1094 | 1280 | −222..+267 | −116..+92 | 0.013 | **0.5536** |
| 01 | N11a (network, all negatives) | 1073 | 1268 | −277..+230 | −187..+82 | 0.017 | 0.5517 |
| 01 | F0 (v2, voxel + 11 features) | 697 | 930 | | | 0.011 | 0.5250 |
| 02 | BASE | 1288 | 1376 | | | 0.008 | 0.5348 |
| 02 | L11 | 1261 | 1380 | −111..+83 | −101..+187 | 0.008 | 0.5306 |
| 02 | N11 | 1100 | 1222 | −300..−64 | −287..−12 | 0.016 | **0.5552** |
| 02 | N11a | 743 | 848 | −914..−403 | −662..−446 | 0.008 | 0.5540 |
| 02 | F0 | 698 | 1037 | | | 0.013 | 0.5271 |

## Reading by the plan

- L11 − N11 at N@1 %: 173 (13.7 % of L11) on raw 01 and 161 (12.8 %) on raw 02. Both are > 5 %, so **the procedure is a bottleneck**.
- N11a does not close the gap (raw 02: 743). Negative subsampling is not the cause.

## Observations (not proven causes)

- The network, with the same 11 inputs as the LDA, ranks entrapment above decoys more (AUC 0.554 / 0.555 against 0.531). Possible cause: the round positives (target winners at q ≤ 0.15, then 0.05) contain false targets. A flexible model with a BCE loss fits these labels and learns any difference between target and decoy library entries. Entrapment entries share that difference with targets. A linear LDA has less freedom to do this.
- F0 (voxel + the same 11 features) is much lower than N11 (697 against 1094). So the voxel branch adds loss on top of the procedure loss.
- L11 (no CNN) is at or above the CNN + LDA base on this population. This needs its own predeclared test before any claim.
