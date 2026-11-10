# Voxel v3 step 0 result

Plan: `PLAN_STEP0.md`, sha256 `916f433e…`, written 2026-11-10T13:40:34Z, before the checks. The sha was checked after the run (unchanged). Output: `STEP0.txt`, exit 0.

## Check A (brightness dependence): NOT confirmed by the frozen rule

| Raw | Score | AUC anchor > random decoys | AUC anchor > brightness-matched decoys | Drop |
|---|---|---|---|---|
| 01 | BASE | 1.0000 | 0.9999 | +0.0001 |
| 01 | V0 | 0.9999 | 0.9987 | +0.0011 |
| 01 | FS | 1.0000 | 0.9999 | +0.0001 |
| 01 | brightness alone | 0.9422 | 0.5030 | +0.4392 |
| 02 | BASE | 1.0000 | 0.9998 | +0.0002 |
| 02 | V0 | 0.9999 | 0.9991 | +0.0007 |
| 02 | FS | 1.0000 | 0.9998 | +0.0002 |
| 02 | brightness alone | 0.9246 | 0.5011 | +0.4236 |

- The matching works: with matched decoys, brightness alone falls to 0.50.
- The V0 drop (0.0011 / 0.0007) is below the 0.01 limit. Non-anchor base hits against matched decoys: V0 0.996 / 0.998.
- Limit of the check: anchors are too easy (all scores give an AUC near 1), so this check cannot see a tail. The ID loss in v2 comes from about 2 in 1000 decoy winners that reach the top. An AUC does not measure a tail of this size. The post-hoc diagnosis (`DIAG_V2f.txt`) showed bright windows enriched about 20× in the decoy top 0.1 %. That is a tail effect and is not a frozen test.

## Check B (voxel score as an extra base LDA feature): FAILS

| Raw | Arm | N@1 % | N@2 % | Δ@1 % vs B0 | Δ@2 % vs B0 | FDP_c q ≤ 0.01 | Entrapment AUC |
|---|---|---|---|---|---|---|---|
| 01 | B0 (12 base features) | 1181 | 1385 | | | 0.013 | 0.5265 |
| 01 | B1 (+ V0) | 976 | 1095 | −358..−41 | −399..−135 | 0.002 | 0.5231 |
| 01 | B2 (+ FS) | 1120 | 1366 | −87..+83 | −100..+56 | 0.013 | 0.5250 |
| 02 | B0 | 1345 | 1463 | | | 0.005 | 0.5293 |
| 02 | B1 (+ V0) | 1217 | 1264 | −236..−84 | −317..−107 | 0.004 | 0.5265 |
| 02 | B2 (+ FS) | 1286 | 1468 | −111..+60 | −96..+38 | 0.009 | 0.5284 |

- With V0 the result is **worse**. The LDA gives V0 the largest weight (+6.3 to +6.5), and the weights of the other features go down. With FS the result does not change.
- Reproduction check failed on raw 01: B0 is 1181 against the population BASE 1024 (+15 %, outside ± 5 %). Raw 02: 1345 against 1288 (+4 %). So B1 / B2 are read only against B0, as the plan says. Possible cause (not tested): the Python LDA trains only on the v2 population (pairs with transition support ≥ 3). That gives 62k decoy negatives instead of 177k in the C++ run.

## Use of the result (by the plan)

A is not confirmed and B fails, so this goes to the user before any v3 training. The brightness cause is not supported by the frozen test. Only the post-hoc tail evidence supports it.

Separate finding, not part of the rule: B0 on the restricted population is above the C++ base on both raws (1181 / 1345 against 1024 / 1288). This needs its own predeclared test, with no change to the gate.
