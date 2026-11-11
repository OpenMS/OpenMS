# Voxel v5 result: fragment mass accuracy (P) and consensus co-elution (C) on top of D11FW

Plan: `PLAN_V5.md` (sha256 `ec7c356dd7…`, written 2026-11-11T13:15:29Z, sha check OK after the run). Base: D11FW (PLAN_V4B). Raws 01 and 02. Entrapment labels are used for evaluation only. DIA-NN is a reference only. Training took about 1.5 min per arm and raw on the laptop (10 threads).

## Verdict: no arm passes; D11FW stays the base

| Arm | Rule 1 (PLAN_V2 1–4 against BASE) | Rule 2 (N@2 % above D11FW, both raws) | Rule 3 (null-winner AUC ≤ D11FW + 0.005) | Result |
|---|---|---|---|---|
| D11FPW (F + P) | PASS: ΔN@2 % +253..+487 / +81..+412; ΔN@1 % −63..+605 / −277..+312; FDP_c 0.013 / 0.014; entrapment AUC 0.492 / 0.496 | FAIL on raw 01 (1657 against 1714); raw 02 1703 against 1620 | PASS (0.4936 / 0.4937) | FAIL |
| D11FPCW (F + P + C) | FAIL: FDP_c at q ≤ 0.01 on raw 01 is 0.017 | FAIL on raw 01 (1664 against 1714) | PASS | FAIL |

## Numbers (decoy-count estimates)

| Raw | Arm | N@0.5 % | N@1 % | N@1.5 % | N@2 % | q ≤ 0.01 N (FDP_c) | Entrapment AUC |
|---|---|---|---|---|---|---|---|
| 01 | D11 | 1256 | 1348 | 1496 | 1572 | 1504 (0.015) | 0.511 |
| 01 | D11FW | 924 | 1521 | 1637 | 1714 | 1628 (0.014) | 0.491 |
| 01 | D11FPW | 1014 | 1478 | 1606 | 1657 | 1538 (0.013) | 0.492 |
| 01 | D11FPCW | 1276 | 1341 | 1492 | 1664 | 1536 (0.017) | 0.495 |
| 02 | D11 | 980 | 1334 | 1417 | 1572 | 1353 (0.012) | 0.514 |
| 02 | D11FW | 885 | 1115 | 1524 | 1620 | 1481 (0.015) | 0.496 |
| 02 | D11FPW | 1013 | 1432 | 1543 | 1703 | 1529 (0.014) | 0.496 |
| 02 | D11FPCW | 830 | 1111 | 1551 | 1684 | 1450 (0.014) | 0.493 |

Paired bootstrap (90 %):

| Raw | Comparison | ΔN@1 % | ΔN@2 % |
|---|---|---|---|
| 01 | D11FPW − D11FW | −43 (−152..+463) | −57 (−104..+47) |
| 02 | D11FPW − D11FW | +317 (−93..+565) | +83 (−128..+159) |
| 01 | D11FPCW − D11FW | −180 (−272..+446) | −50 (−213..+49) |
| 02 | D11FPCW − D11FW | −4 (−314..+485) | +64 (−95..+201) |

Funnel at 1 % (DIA-NN reference):

| Raw | Arm | C wrong window | D right window, below cut | E accepted |
|---|---|---|---|---|
| 01 | D11FW | 113 | 847 | 1432 |
| 01 | D11FPW | 88 | 909 | 1395 |
| 02 | D11FW | 109 | 1213 | 1078 |
| 02 | D11FPW | 87 | 932 | 1381 |

## Observations (not part of the verdict)

1. **Mass accuracy steadies the 1 % cut on raw 02.**
   - D11FPW N@1 % on raw 02 is 1432, against 1115 for D11FW and 1334 for D11. Stage E rises by 303 there.
   - Raw 01 loses 37 at 1 % and 57 at 2 %.
   - Across the two raws, D11FPW gives the most balanced curve of any arm so far: N@1 % 1478 / 1432, N@2 % 1657 / 1703.
   - The paired intervals against D11FW include 0 everywhere.
2. **Consensus co-elution (C) does not help.** It shifts the 1 % cut on raw 01 down (N@1 % 1341) and raises FDP_c to 0.017.
3. **Rule 2 is strict against D11FW.** D11FW has the highest raw 01 N@2 % so far (1714), and its raw 02 1 % weakness is not tested by rule 2. A rule on both 1 % and 2 % (or the paired-bootstrap rule pending user approval) would rank D11FPW differently. Changing the rule needs user approval.
