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
