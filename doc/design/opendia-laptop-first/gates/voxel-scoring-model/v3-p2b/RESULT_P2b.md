# Voxel v3 phase 2b result

Plan: `PLAN_V3_P2b.md` (sha `f40e4cb9…`, written 2026-11-10T15:28:25Z, before any run, checked after gate 1).

## Gate 1 (11 features, one network, decoy-whitened objective): PASS

| Raw | Arm | N@0.5/1/1.5/2 % | FDP_c q ≤ 0.01 | q ≤ 0.01 N | Entrapment AUC | Gate 1 (N@1 % ≥ 0.95 × L11d, AUC ≤ L11d + 0.005) |
|---|---|---|---|---|---|---|
| 01 | BASE (CNN + LDA) | 900 / 1024 / 1267 / 1314 | 0.013 | 1081 | 0.5313 | |
| 01 | L11d (closed form, reference) | 1036 / 1259 / 1365 / 1496 | 0.013 | 1285 | 0.5274 | |
| 01 | **D11** | **1256 / 1348 / 1496 / 1572** | 0.015 | 1504 | **0.5109** | PASS (1348 ≥ 1196; 0.5109 ≤ 0.5324) |
| 01 | D11n (no z skip, information only) | 525 / 528 / 577 / 580 | 0.004 | 510 | 0.5143 | |
| 02 | BASE | 1215 / 1288 / 1333 / 1376 | 0.008 | 1241 | 0.5348 | |
| 02 | L11d | 1116 / 1244 / 1299 / 1432 | 0.007 | 1215 | 0.5273 | |
| 02 | **D11** | **980 / 1334 / 1417 / 1572** | 0.012 | 1353 | **0.5141** | PASS (1334 ≥ 1182; 0.5141 ≤ 0.5323) |
| 02 | D11n | 318 / 381 / 1400 / 1435 | 0.006 | 327 | 0.5131 | |

- D11 against the population BASE: ΔN@2 % 90 % interval +165..+400 (raw 01) and +16..+305 (raw 02), both above 0. ΔN@1 %: +81..+519 and −344..+218.
- The entrapment AUC is the lowest so far (0.511 / 0.514, base 0.531 / 0.535). So the learned features separate "present" from "absent", not "target entry" from "decoy entry".
- FDP_c at q ≤ 0.01 is 0.015 on raw 01, exactly on the PLAN_V2 limit. The decoy-count estimate is slightly optimistic there.
- D11n (no skip of the raw features) is unstable (raw 02 N@1 % 381 but N@1.5 % 1400). The skip is needed.
- These are heuristic-only results. D11 is not a voxel model and is not part of the gate 2 verdict.

## Gate 2 (voxel encoder + 11 features, one network, same objective): FAIL

Plan sha checked after gate 2 (OK). Output: `EVALUATION_G2.txt`. Training took 29 min (raw 01) and 22 min (raw 02), with 4.1 GB peak memory.

| Raw | Arm | N@0.5/1/1.5/2 % | Δ@1 % vs BASE | Δ@2 % vs BASE | FDP_c q ≤ 0.01 | q ≤ 0.01 N | Entrapment AUC |
|---|---|---|---|---|---|---|---|
| 01 | BASE | 900 / 1024 / 1267 / 1314 | | | 0.013 | 1081 | 0.5313 |
| 01 | D11 | 1256 / 1348 / 1496 / 1572 | +73..+505 | +158..+396 | 0.015 | 1504 | 0.5109 |
| 01 | **DV** | 1254 / 1388 / 1447 / 1517 | +66..+521 | +121..+380 | **0.019** | 1497 | **0.4996** |
| 02 | BASE | 1215 / 1288 / 1333 / 1376 | | | 0.008 | 1241 | 0.5348 |
| 02 | D11 | 980 / 1334 / 1417 / 1572 | −348..+224 | +16..+311 | 0.012 | 1353 | 0.5141 |
| 02 | **DV** | 1027 / 1291 / 1394 / 1438 | −147..+142 | **−51..+282** | 0.005 | 1131 | **0.4984** |

(The D11 intervals differ slightly from gate 1 because the bootstrap stream depends on the arm order. The N values are the same.)

| Rule (PLAN_V2, both raws) | Raw 01 | Raw 02 |
|---|---|---|
| 1. Δ@2 % interval above 0 | pass | **fail** (−51..+282) |
| 2. Δ@1 % interval not entirely below 0 | pass | pass |
| 3. FDP_c at q ≤ 0.01 ≤ 0.015 | **fail** (0.019) | pass |
| 4. Entrapment AUC ≤ BASE | pass | pass |

**DV fails** (rule 1 on raw 02, rule 3 on raw 01).

## Observations (not proven causes)

- The objective fixes the v2 failure mode. The voxel model is no longer below the base: DV gives 1388 / 1291 at 1 % (v2 FS: 921 / 1040), and its entrapment AUC is 0.500 / 0.498. This is the first model that does not rank entrapment above decoys.
- The voxel input adds no IDs over the 11 features alone. DV − D11 at N@2 %: −55 (raw 01) and −134 (raw 02). At N@1 %: +40 and −43.
- Raw 01 DV FDP_c at q ≤ 0.01 is 0.019. The decoy-count q is optimistic there, even though the entrapment AUC is 0.50. Possible cause: the top of the score is tight, and a few entrapment hits fall just above the decoys. Not tested.
- Without the voxel input, D11 meets all four PLAN_V2 numbers on both raws: Δ@2 % +158..+396 and +16..+311, FDP_c 0.015 / 0.012, AUC 0.511 / 0.514. But D11 was gate 1 of this plan. The plan did not name it as a candidate for the voxel-model gate, so this is an observation, not a pass. It is also not a voxel model.
