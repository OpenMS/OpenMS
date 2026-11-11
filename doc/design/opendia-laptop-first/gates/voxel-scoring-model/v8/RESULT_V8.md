# Voxel v8 result: D11FWT on the linear-end-RTMap population. FAIL (rule 3 on both raws)

**Plan.** `PLAN_V8.md`, sha256 `cf14d53d3f…`, written 2026-11-11T15:56:15Z, before any rtlin score was computed. The sha was checked OK by `chain_score.sh` before training.

**Data.** `rtlin-20261111/chain_data.sh`, run with the binaries in `binary.sha256` (commit `6015a6e53f`).
- The 6-fragment re-export reproduced the search features byte for byte on both raws (`CMP6.txt` IDENTICAL).
- The 12-fragment XIC check passes: median relative difference 0, p99 0.031.
- **Deviation:** the laptop suspended from 17:09 to 18:11 UTC during `export6_02`. The processes resumed, nothing was restarted and the outputs were checked as above. See `NOTES.txt`.

## Calibration (same options as before)

| Raw | Knot library-RT range, old → new | Anchors used | RT sigma (s) |
|---|---|---|---|
| 01 | 0.098–0.752 → 0.066–0.796 | 504 → 527 | 21.05 → 21.02 |
| 02 | 0.098–0.773 → 0.051–0.796 | 596 → 554 | 20.74 → 19.54 |

## Population

| Raw | Windows, old → new | Pairs, old → new |
|---|---|---|
| 01 | 1,270,510 → 1,155,776 | 94,109 → 85,422 |
| 02 | 1,275,383 → 1,038,687 | 95,447 → 81,329 |

## Result (`EVALUATION_V8.txt`, `NULL_V8.txt`)

| Raw | Arm | N@FDP_c 0.5 / 1 / 1.5 / 2 % | q ≤ 0.01 N, FDP_c | AUC (entrapment > decoy) | 90 % vs BASE @1 % | 90 % vs BASE @2 % |
|---|---|---|---|---|---|---|
| 01 | BASE (rebuilt) | 873 / 1105 / 1293 / 1337 | 1102, 0.009 | 0.5352 | | |
| 01 | D11FWT | 1142 / **1519** / 1694 / **1839** | 1728, **0.016** | 0.4919 | +92..+693 | +336..+591 |
| 02 | BASE (rebuilt) | 939 / 1033 / 1168 / 1205 | 866, 0.002 | 0.5295 | | |
| 02 | D11FWT | 1250 / **1420** / 1644 / **1736** | 1716, **0.019** | 0.4992 | +125..+638 | +337..+674 |

**Rules on both raws:**
- Rules 1, 2 and 4 pass.
- **Rule 3 (FDP_c at q ≤ 0.01 ≤ 0.015) fails:** 0.016 on raw 01 and 0.019 on raw 02.
- Null-winner AUC is 0.493 / 0.495, against BASE 0.522 / 0.521.

## Reported only

**Old population D11FWT:** N@1 % 1496 / 1426, N@2 % 1716 / 1736, q ≤ 0.01 FDP_c 0.0137 / 0.0131. New: N@1 % 1519 / 1420, N@2 % 1839 / 1736. The populations differ, so this is not a paired comparison.

**Funnel against DIA-NN** (`funnel_D11FWT.txt`), raw 01 / raw 02:

| Stage | Old | New |
|---|---|---|
| A | 185 / 151 | 133 / 120 |
| B | 166 / 180 | 97 / 83 |
| C | 113 / 109 | 119 / 112 |
| D | 866 / 929 | 954 / 1060 |
| E | 1413 / 1362 | 1440 / 1356 |

A+B falls by 121 / 128. Most of the recovered precursors land in D, not in E.

**The q ≤ 0.01 cut is decoy-calibrated but not entrapment-calibrated** (descriptive):

| Raw | Old: decoys / entrapment | New: decoys / entrapment | Targets, new |
|---|---|---|---|
| 01 | 15 / 11 | 16 / 14 | 1728 |
| 02 | 14 / 10 | 16 / 16 | 1716 |

- If false targets behaved like decoys, about 7.9 entrapment would be expected in each new set.
- The chance of the observed count or more under that expectation is 0.03 (raw 01) and 0.008 (raw 02).
- The binomial 90 % interval of FDP_c is 0.010–0.026 (raw 01) and 0.012–0.029 (raw 02).
- The counts are small. Their direction matches V7 (target/decoy asymmetry), not noise alone.

**Where the excess sits** (`DIAG_V8_ENDS.txt`): by library RT, inside the old knot range against the new ends, at q ≤ 0.01:
- **Ends:** 202 / 207 targets, 2 / 1 decoys, 1 / 2 entrapment. The ends add about 80–100 accepted targets per raw at a low decoy count.
- **Inside the old knot range:** 1526 / 1509 targets, 14 / 15 decoys, 13 / 14 entrapment. The entrapment-to-decoy ratio is 0.93 against 0.71 in the old population.
- So the higher FDP_c is **not** shown to come from the extrapolated ends. It sits in the unchanged part of the RT range. A plausible reading, not shown, is that the new rebuild changes the training population (seeds, anchors, BASE features), and the model's decoy-relative calibration there is less entrapment-faithful.

## Outcome (per PLAN_V8)

- **Fail.** D11FWT does not become the base on the rtlin population.
- The C++ RTMap change is **not** reverted. It is kept as committed, pending the user's decision.
- Open questions for the user:
  1. Keep the linear-end RTMap. It recovers A+B and its end region is decoy-clean.
  2. Look for the cause of the inside-range entrapment excess. The decoy and entrapment counts at q ≤ 0.01 are only 10–16.
