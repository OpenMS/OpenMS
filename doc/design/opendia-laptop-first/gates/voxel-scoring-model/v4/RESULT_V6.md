# Voxel v6 result: per-fold decoy-tail calibration of the final score

Plan: `PLAN_V6.md` (sha256 `73258acfb9…`, written 2026-11-11T13:44:16Z, sha check OK after the run). Base: D11FW. Raws 01 and 02. There is no retraining: `calib_tail.py` is a label-free, monotone, per-fold transform built only from decoy winners. Checks in `calib_tail.log`: order within each fold and window choice are unchanged. Entrapment labels are used for evaluation only.

## Verdict: D11FWT passes; D11FPWT fails

| Arm | Rule 1 (PLAN_V2 1–4 against BASE) | Rule 2 (N@2 % above D11FW, both raws) | Rule 3 (null-winner AUC ≤ D11FW + 0.005) | Result |
|---|---|---|---|---|
| D11FWT (D11FW + tail calibration) | PASS: ΔN@2 % +293..+505 / +162..+448; ΔN@1 % −176..+659 / −99..+380; FDP_c 0.0137 / 0.0131; entrapment AUC 0.491 / 0.496 | PASS (1716 > 1714; 1736 > 1620) | PASS (0.4933 / 0.4936) | **PASS** |
| D11FPWT (D11FPW + tail calibration) | FAIL: FDP_c at q ≤ 0.01 on raw 01 is 0.0164 | FAIL on raw 01 (1677 < 1714) | PASS | FAIL |

**Rule 2 is passed narrowly on raw 01:** 1716 against 1714, with a paired interval of −22..+27.

## Numbers (decoy-count estimates; FDP_c at q ≤ 0.01 exact)

| Raw | Arm | N@0.5 % | N@1 % | N@1.5 % | N@2 % | q ≤ 0.01 N (FDP_c) | Entrapment AUC | Null-winner AUC |
|---|---|---|---|---|---|---|---|---|
| 01 | BASE | 900 | 1024 | 1267 | 1314 | 1081 (0.013) | 0.531 | 0.518 |
| 01 | D11 | 1256 | 1348 | 1496 | 1572 | 1504 (0.0147) | 0.511 | 0.5015 |
| 01 | D11FW | 924 | 1521 | 1637 | 1714 | 1628 (0.0136) | 0.491 | 0.4936 |
| 01 | **D11FWT** | 971 | 1496 | 1616 | 1716 | 1614 (0.0137) | 0.491 | 0.4933 |
| 01 | D11FPWT | 1275 | 1490 | 1571 | 1677 | 1600 (0.0164) | 0.492 | 0.4934 |
| 02 | BASE | 1215 | 1288 | 1333 | 1376 | 1241 (0.008) | 0.535 | 0.523 |
| 02 | D11 | 980 | 1334 | 1417 | 1572 | 1353 (0.0119) | 0.514 | 0.5022 |
| 02 | D11FW | 885 | 1115 | 1524 | 1620 | 1481 (0.0150) | 0.496 | 0.4934 |
| 02 | **D11FWT** | 1127 | 1426 | 1625 | 1736 | 1536 (0.0131) | 0.496 | 0.4936 |
| 02 | D11FPWT | 1325 | 1459 | 1635 | 1687 | 1519 (0.0119) | 0.496 | 0.4937 |

Paired bootstrap against D11FW (90 %; `VS_D11FW_V6.txt`):

| Raw | Arm | ΔN@1 % | ΔN@2 % |
|---|---|---|---|
| 01 | D11FWT | −25 (−28..+153) | +2 (−22..+27) |
| 02 | D11FWT | **+311 (+2..+696)** | +116 (−9..+207) |
| 01 | D11FPWT | −31 (−127..+517) | −37 (−105..+69) |
| 02 | D11FPWT | +344 (−28..+696) | +67 (−64..+238) |

Funnel at 1 % (DIA-NN reference):

| Raw | Arm | C wrong window | D right window, below cut | E accepted |
|---|---|---|---|---|
| 01 | D11FW | 113 | 847 | 1432 |
| 01 | D11FWT | 113 | 866 | 1413 |
| 02 | D11FW | 109 | 1213 | 1078 |
| 02 | D11FWT | 109 | 929 | 1362 |

## Reading

- **What the calibration fixes.** It removes the raw 02 1 % instability of D11FW. N@1 % rises from 1115 to 1426, and DIA-NN-confirmed precursors at 1 % rise from 1078 to 1362. FDP_c at q ≤ 0.01 falls from 0.0150 to 0.0131. On raw 01 nothing changes within noise.
- **Why it likely works.** The probable cause (consistent with the pre-plan look, not proven) was the pooling of folds whose decoy tails differ, after a median/MAD scaling that fits only the bulk.
  - The calibration maps every fold's score to the same decoy-tail scale.
  - It uses no labels and keeps each fold's order.
  - It is also what a single model in C++ would need when it scores folds (or runs) on one scale.
- **D11FW against D11FWT.** D11FWT is the new base under PLAN_V6. Against D11 (PLAN_V3_P2b), D11FWT is +148 / +92 at 1 % and +144 / +164 at 2 %.
- **Limits.**
  - Only two raws; raw 03 is not used.
  - The raw 01 rule-2 margin is +2.
  - Mass accuracy (P) with calibration fails on raw 01 FDP_c (13 entrapment winners at q ≤ 0.01).
