# Voxel v2 result: FAIL, no arm passes

Plan: `PLAN_V2.md`, sha256 `d156805d6b9e724ca62411be30df35bf17dd7d8a130dcd0d571780e2e62cb787` (checked after the run, unchanged), written 2026-11-10T10:52:28Z, before training.
Chain: `chain_v2.sh`, finished 2026-11-10T12:22Z. All steps exit 0.

## Checks before training

- Exports raw 01 and raw 02: `search.features.tsv` is byte-identical to the base on both raws (`CMP.txt` = IDENTICAL, `EXIT` = 0).
- Raw 01 export precursor list = `raw01/export_precursors.txt` (sorted diff: only the header comment is different).
- Voxels: all windows have fragment traces, 98 % of windows have MS1 points, 99.997 % of windows match the library.

## Result (full output: `EVALUATION_V2.txt`)

N = targets at the given FDP_c. The intervals are bootstrap 90 % intervals of the arm minus the population base.

| Raw | Arm | N@1 % | N@2 % | Δ@1 % interval | Δ@2 % interval | FDP_c at q ≤ 0.01 | Entrapment AUC |
|---|---|---|---|---|---|---|---|
| 01 | BASE | 1024 | 1314 | | | 0.013 | 0.5313 |
| 01 | V0 | 422 | 524 | −1097..−386 | −965..−545 | 0.000 | 0.5117 |
| 01 | F0 | 697 | 930 | −640..−83 | −538..−226 | 0.011 | 0.5250 |
| 01 | S0 | 638 | 706 | −963..−217 | −729..−343 | 0.008 | 0.5202 |
| 01 | FS | 921 | 1179 | −338..+140 | −267..+17 | **0.018** | **0.5354** |
| 02 | BASE | 1288 | 1376 | | | 0.008 | 0.5348 |
| 02 | V0 | 278 | 718 | −1226..−500 | −1056..−491 | 0.011 | 0.5182 |
| 02 | F0 | 698 | 1037 | −678..−252 | −638..−224 | 0.013 | 0.5271 |
| 02 | S0 | 519 | 768 | −1131..−486 | −901..−431 | 0.006 | 0.5080 |
| 02 | FS | 1040 | 1235 | −477..−42 | −288..−30 | 0.007 | **0.5366** |

## Rule check

| Arm | Rule 1 (Δ@2 % > 0) | Rule 2 (Δ@1 % not all < 0) | Rule 3 (FDP_c ≤ 0.015) | Rule 4 (AUC ≤ base) | Verdict |
|---|---|---|---|---|---|
| V0 | fail 01, 02 | fail 01, 02 | pass | pass | FAIL |
| F0 | fail 01, 02 | fail 01, 02 | pass | pass | FAIL |
| S0 | fail 01, 02 | fail 01, 02 | pass | pass | FAIL |
| FS | fail 01, 02 | pass 01, fail 02 | fail 01 (0.018) | fail 01, 02 | FAIL |

## Observations (not proven causes)

- The voxel model alone (V0) is much weaker than the base on the full population: 422 / 278 at 1 % against 1024 / 1288. The pilot sample did not show this gap.
- The heuristic features give the largest gain (V0 to F0). The rounds add more on top of the features (F0 to FS), but the rounds alone help little (V0 to S0).
- FS is the only arm near the base. But with FS, the entrapment AUC is above the base on both raws, and FDP_c at q ≤ 0.05 goes up (0.079 / 0.085 against 0.061 / 0.074). A possible cause: the rounds take target winners as positives, and some of these are entrapment hits, so the model can learn to prefer them.
- Decoy-like windows are negatives only in round 0. The rounds use decoy winners as negatives.

## Next step

No arm passes, so by the plan the user decides the next step. The direction that was agreed before is v3: a lower weight for the decoy-like negatives. A note for this decision: FS already prefers entrapment more than the base does. A lower weight on decoy-like negatives can make this worse, because these negatives are the windows most like an absent target.
