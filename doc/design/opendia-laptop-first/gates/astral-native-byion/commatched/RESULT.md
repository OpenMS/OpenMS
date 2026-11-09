# k11 result: commonness-matched decoys, native OpenDIA, raws 01 and 02. FAIL under the PLAN.md rule.

| run | library | N at FDP_c 0.5/1/1.5/2% | q<=0.01 N / FDP_c | q<=0.05 N / FDP_c |
|---|---|---|---|---|
| raw 01 | base | 899 / 1029 / 1270 / 1321 | 1068 / 0.013 | 1638 / 0.065 |
| raw 01 | commonness-matched | 1023 / 1101 / 1219 / 1334 | 1033 / 0.006 | 1663 / 0.061 |
| raw 02 | base | 1230 / 1298 / 1331 / 1372 | 1236 / 0.007 | 1788 / 0.074 |
| raw 02 | commonness-matched | 1178 / 1210 / 1293 / 1316 | 1213 / 0.002 | 1672 / 0.060 |

- Library (build.txt): 26% of the chosen decoys are pseudo-reverse and 74% are shuffles. Decoy commonness mean goes from 0.668 to 0.730 (targets 0.792). The AUC of target against decoy commonness goes from 0.580 to 0.536.
- q<=0.05 FDP_c: raw 01 goes from 0.065 to 0.061 (drop 0.005). Raw 02 goes from 0.074 to 0.060 (drop 0.014). Both values are above 0.06 before rounding (EVALUATION.txt, condition False).
- q<=0.01 FDP_c is more conservative on both raws (0.006 and 0.002).
- N at FDP_c 1%: raw 01 increases by 7% (1101). Raw 02 decreases by 7% (1210), below the 95% limit.
- Reading: FAIL. The 5% drop is 0.01 or more on raw 02 only, and raw 02 N falls below 95% of base. The decoys move the tail in the correct direction but do not close the gap. One limit: the chosen decoys still have a commonness mean 0.06 below the targets.
- Cost: crossfit 28 min (raw 01) and 22 min (raw 02), 25 GB peak. Search 3.4 min and 5.7 min.
