# RESULT: effect of the decoy library on OpenDIA and DIA-NN IDs (Astral entrapment raws 01/02, PAIRLIB targets)

Targets are identical in every library: 302,480 precursors, half of them entrapment. Metrics, all evaluation-only:
- **q ≤ 0.01 N:** the IDs the tool reports at its own q ≤ 0.01.
- **FDP_c:** the entrapment-estimated false discovery proportion among them.
- **N@FDP_c 1 %:** IDs at a true (entrapment) FDP of 1 %. This is a ranking measure that does not depend on how the
  tool's q-values are calibrated.

## OpenDIA (D11FDWT; codon runs of PAIRLIB and PLAN_DECOYGEN)
| decoys | q ≤ 0.01 N, raw01 / raw02 | FDP_c | N@FDP_c 1 % | N@FDP_c 2 % |
|---|---|---|---|---|
| old: pseudo-rev + K/R swap, copied assay | 1978 / 1784 | 0.032 / 0.017 | 1549 / 1671 | 1821 / 1870 |
| D: same sequences, own predicted assay (DIALibGen) | 1667 / 1445 | 0.008 / 0.008 | 1733 / 1618 | 1864 / 1801 |
| E: D built by C++ `OpenSwathDecoyGenerator -predict_decoys` | 1667 / 1457 | 0.010 / 0.007 | 1694 / 1628 | 1908 / 1903 |
| A: shuffle, own predicted assay | 837 / 1545 | 0.002 / 0.010 | 1639 / 1678 | 1917 / 1880 |
| B: DIA-NN-style mutate, own predicted assay | 1425 / 1182 | 0.004 / 0.003 | 1716 / 1653 | 1860 / 1868 |
| C: pseudo-rev without swap, own predicted assay | 1689 / 1617 | 0.011 / 0.009 | 1654 / 1659 | 1949 / 1846 |

## DIA-NN 2.2 (same settings as diann-pairlib-20261113; library decoys used: 604,480 / 598,279 library precursors loaded)
| decoys | q ≤ 0.01 N, raw01 / raw02 | FDP_c | q ≤ 0.05 FDP_c | N@FDP_c 1 % | N@FDP_c 2 % |
|---|---|---|---|---|---|
| DIA-NN's own (targets only) | 2361 / 2322 | 0.014 / 0.010 | 0.060 / 0.060 | 2067 / 2301 | 2522 / 2543 |
| old library decoys | 2064 / 2027 | 0.009 / 0.008 | 0.033 / 0.043 | 2154 / 2108 | 2354 / 2310 |
| D library decoys (own predicted assay) | 1909 / 1656 | 0.007 / 0.005 | 0.037 / 0.032 | 2090 / 1970 | 2203 / 2184 |

The raw01 own-decoy row is the 2026-11-13 run, reproduced exactly by `evaluate.py`.

Runtime and memory:
- Library-decoy runs: 20–21 min each with 9 threads, two in parallel.
- Own-decoy run on raw02: 8:48.
- The first attempt with `--report-decoys` crashed DIA-NN 2.2 (exit 139) while it wrote the report, so the runs were
  repeated without that option.

## Reading
- **OpenDIA:**
  - The predicted decoys fix the calibration: FDP_c at q ≤ 0.01 goes from 0.032 / 0.017 to 0.008 / 0.008.
  - Reported IDs at q ≤ 0.01 fall by about 310–340 per raw, because the old q-values were optimistic.
  - Ranking quality (N@FDP_c 1 %) moves +184 on raw01 and −53 on raw02. Its noise between retrainings is not measured,
    so these moves are not a clear gain or loss. The C++ arm E is within the same range.
- **DIA-NN:**
  - It was already calibrated, or conservative, with the old library decoys (FDP_c 0.009 / 0.008). So the
    under-counting with old decoys is specific to OpenDIA's scorer, not to the decoys alone.
  - With the predicted decoys it becomes more conservative: FDP_c 0.007 / 0.005 at q ≤ 0.01, and about 0.03 at
    q ≤ 0.05. It reports fewer IDs (−155 / −371 against the old decoys).
  - Its ranking is slightly lower (N@FDP_c 1 %: −64 / −138).
  - Its own decoys still give the most IDs (2361 / 2322), at FDP_c 0.014 / 0.010, slightly above nominal at q ≤ 0.01
    and about 1.2× at q ≤ 0.05.
- **Gap at matched true FDP:** at N@FDP_c 1 %, OpenDIA with D or E reaches 1618–1733. DIA-NN reaches 2067 / 2301 with
  its own decoys and 1970–2154 with any library decoys. The ranking gap of about 350–680 IDs is unchanged by the decoys.
