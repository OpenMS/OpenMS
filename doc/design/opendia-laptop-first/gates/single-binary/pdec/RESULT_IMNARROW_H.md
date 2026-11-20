# RESULT_IMNARROW_H: narrower 1/K0 window, held-out run H (reported, 20 Nov 2026)

H (22474.d) was submitted at 08:58:49 UTC, after the D verdict was committed (ab765855cc, 08:58:46). Per
PLAN_IMNARROW it is reported, not gated. The setup is the same as for D: pd10 0942ffbdf3, the PLAN_ODP predicted
library, CUR = default, NEW = `-im_tolerance_scale 0.6`, 5 seeds. Output is in `EVAL_IMNARROW_H.txt` and
`COST_IMNARROW_H.txt`.

| seed | CUR N at q ≤ 0.01 / E / N at FDP_c 1 % | NEW N at q ≤ 0.01 / E / N at FDP_c 1 % |
|---|---|---|
| 0 | 26,224 / 281 / 25,594 | 27,274 / 300 / 26,697 |
| 1 | 26,317 / 288 / 25,657 | 27,280 / 289 / 26,577 |
| 2 | 26,357 / 289 / 25,578 | 27,308 / 292 / 26,647 |
| 3 | 26,311 / 279 / 25,759 | 27,071 / 273 / 26,651 |
| 4 | 26,246 / 276 / 25,643 | 27,180 / 287 / 26,522 |

| | CUR | NEW |
|---|---|---|
| Median N at FDP_c 1 % | 25,643 | **26,647 (+3.9 %)** |
| Pooled lower bound at q ≤ 0.01 (90 %) | 1.07 % (1.03-1.12) | 1.06 % (1.01-1.11) |
| Pooled combined bound | 1.23 % | 1.21 % |
| Pair win | 0.483 | 0.492 |
| Median wall / peak RSS | 11:10 / 35.1 GB | 10:51 / 33.5 GB |

- H reproduces D: +3.9 % here against +4.3 % on D.
- Every NEW seed beats every CUR seed: the lowest NEW seed has 26,522 and the highest CUR seed 25,759.
- Calibration is unchanged, and there is no extra cost.
- The default is now 0.6 (OpenMS 9ce3370ee3). That default is byte-identical to the tested NEW arm on r ≈ 1.56 HeLa D s0, `-im_tolerance_scale 1` is byte-identical to the old default, and 25/25 ctests pass.
- Against DIA-NN on this library (26,006 on D, 41,272 on H), OpenDIA is now at 45 % on D (11,815) and 65 % on H (26,647).
