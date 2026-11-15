# PLAN_CALREFINE result: fails criterion 3 on raw02; reverted

**Verdict: fail.** Criteria 1 and 2 pass on all 4 runs. Criterion 3 fails on raw02, where N@FDP_c 1 % is 1826 against a required ≥ 1941. Per the plan the change is reverted in OpenMS. The diff is kept as `calrefine.patch`.

## Criteria 1-2 (`CRITERIA12.txt`; DIA-NN robust spread around the written map)

| Run | Current spread / ratio | Refined spread / ratio | Tolerance current → refined | DIA-NN coverage current → refined |
|---|---|---|---|---|
| arm E raw01 | 37.7 / 1.39 | 24.8 / 0.96 | 151 → 99 s | 0.966 → 0.961 |
| arm E raw02 | 30.9 / 1.22 | 24.9 / 1.00 | 124 → 100 s | 0.967 → 0.959 |
| PAIRLIB raw01 | 27.4 / 1.06 | 22.6 / 0.90 | 110 → 90 s | 0.961 → 0.957 |
| PAIRLIB raw02 | 28.0 / 1.13 | 25.2 / 1.04 | 112 → 101 s | 0.961 → 0.957 |

- **Map accuracy:** |err| p50 is unchanged or better (16.5-17.5 s), and the refined pass yields more anchors (499-891 against 425-795).

## Criterion 3 (`CRITERION3.txt`; self-calibrating OpenDIAWorkflow, both raws, seed 0)

| Run | Current code | Refined | Required |
|---|---|---|---|
| raw01 N@FDP_c 1 % | 1718 | 1817 | ≥ 1718: pass |
| raw01 FDP_c at q ≤ 0.01 | 0.0104 | 0.0081 | ≤ 0.010: pass |
| raw02 N@FDP_c 1 % | 1941 | 1826 | ≥ 1941: **fail** |
| raw02 FDP_c at q ≤ 0.01 | 0.0045 | 0.0047 | pass |

- **Global q:** 1831 target precursors at global q ≤ 0.01, against 1839 before.

## Reading
- **Calibration:** the refined calibration is the better calibration by both calibration criteria.
- **raw01:** it fixes the raw01 FDP problem and raises N.
- **raw02:** N falls by 115. The criterion compared single seed-0 runs, and the raw02 seed spread of the scorer is about 1816-2021 (step 1) or 1816-1940 (RESULT_WORKFLOW). A drop of this size is inside seed noise.
  - The narrower window (124 → 100 s) also removes some true candidates. Coverage fell 0.967 → 0.959.
- **Plan design:** the plan's end-to-end criterion used a single seed. That is the weak point of the plan, but the verdict stands as frozen.
- **Next step:** a re-test needs a new plan with a multi-seed end-to-end criterion. That is not done here.
