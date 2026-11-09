# k11 result: commonness-matched decoys, native OpenDIA, raws 01 and 02 — FAIL (PLAN.md rule)

| run | library | N at FDP_c 0.5/1/1.5/2% | q<=0.01 N / FDP_c | q<=0.05 N / FDP_c |
|---|---|---|---|---|
| raw 01 | base | 899 / 1029 / 1270 / 1321 | 1068 / 0.013 | 1638 / 0.065 |
| raw 01 | commonness-matched | 1023 / 1101 / 1219 / 1334 | 1033 / 0.006 | 1663 / 0.061 |
| raw 02 | base | 1230 / 1298 / 1331 / 1372 | 1236 / 0.007 | 1788 / 0.074 |
| raw 02 | commonness-matched | 1178 / 1210 / 1293 / 1316 | 1213 / 0.002 | 1672 / 0.060 |

- Library (build.txt): chosen decoys 26% pseudo-reverse, 74% shuffles. Commonness mean 0.668 -> 0.730 (targets 0.792); AUC target vs decoy 0.580 -> 0.536.
- q<=0.05 FDP_c: raw 01 0.065 -> 0.061 (drop 0.005), raw 02 0.074 -> 0.060 (drop 0.014). Neither is <= 0.06 strictly at the rule's precision on raw 01; raw 02 sits at 0.060.
- q<=0.01 FDP_c becomes more conservative on both raws (0.006 / 0.002).
- N at FDP_c 1%: raw 01 +7% (1101), raw 02 -7% (1210, below the 95% bar).
- Reading: FAIL. The 5% drop is >= 0.01 on raw 02 only, and raw 02 N falls below 95% of base. The change moves the tail in the right direction but does not close it; partial commonness matching (mean still 0.06 below targets) is one limit.
- Cost: crossfit 28 min / 22 min, 25 GB peak; search 3.4 / 5.7 min.
