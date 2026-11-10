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

## Phase 1b (PLAN_V3_P1b.md, sha `2dea081f…`, written 14:52:57Z, checked after the run): linear head with BCE (G11)

| Raw | Arm | N@1 % | N@2 % | FDP_c q ≤ 0.01 | Entrapment AUC |
|---|---|---|---|---|---|
| 01 | L11 (LDA) | 1267 | 1484 | 0.010 | 0.5310 |
| 01 | N11 (MLP, BCE) | 1094 | 1280 | 0.013 | 0.5536 |
| 01 | G11 (linear, BCE) | 1300 | 1345 | 0.004 | **0.5577** |
| 02 | L11 | 1261 | 1380 | 0.008 | 0.5306 |
| 02 | N11 | 1100 | 1222 | 0.016 | 0.5552 |
| 02 | G11 | 1146 | 1242 | 0.002 | **0.5571** |

Reading by the plan: the "flexible head" branch fails. Raw 02 G11 is 9.1 % below L11, and on both raws the G11 entrapment AUC is 0.026 above L11 (limit 0.005). Raw 02 matches the "loss / schedule" branch (gap > 5 %). Raw 01 matches neither branch, because G11 ≥ L11 at 1 %.

The clearest signal: a **linear** model with the same 11 inputs prefers entrapment as much as the MLP does when it is trained with BCE on the round labels (0.557 against 0.531 for the LDA). So the entrapment preference comes from the training objective or the sampled negatives, not from model capacity. Possible cause (not proven): BCE puts its weight on the boundary examples. The round positives near the boundary are where the false targets are, so BCE fits them. LDA uses class means and a pooled covariance and is less driven by the boundary. A sampled-negative effect is not excluded (G11 uses 6000 sampled negatives, L11 uses all).

Next (needs its own frozen plan): keep one network, change only the objective. Candidates: a Fisher-discriminant (LDA) loss on the last layer, or a paired target/decoy ranking loss. Test first on the 11 features (it must match L11, AUC included), then add the voxel encoder.
