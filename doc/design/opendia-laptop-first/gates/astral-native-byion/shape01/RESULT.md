# Stage-2 shape group: raw 01 native result: FAIL (calibration), large ID gain

`PLAN.md` was written before the run. Binary is commit `f98c0dca2c`. Columns 1–6 and 9 of crossfit, and feature columns 1–16, are byte-identical to the byion run, so the picker refactor did not change the boundary features (`CMP.txt`).

The metric is `eval_arm.py`: entrapment-matched N ranked by the per-fold decoy-standardised score. Under this metric the byion reference gives 1010; the earlier `evaluate.py` ranked by raw score and gave 1029.

| | FDP_c 0.5% | **1%** | 1.5% | 2% | plain q≤0.01 / 0.05 (FDP_c) | bulk-null q≤0.01 / 0.05 (FDP_c, 3 seeds) |
|---|---|---|---|---|---|---|
| byion (ref) | 925 | 1010 | 1248 | 1335 | 1068 (0.013) / 1638 (0.065) | 939 (0.006) / 1520–1540 (0.048–0.050) |
| **+shape** | 1056 | **1233 (+22%)** | 1267 | 1445 | 1227 (0.010) / 1867 (**0.091**) | 1196–1213 (0.010) / 1717–1737 (**0.068–0.072**) |
| OSW+Percolator | 1015 | 1151 | 1251 | 1286 | | |

- **N criterion: PASS.** +22% at FDP_c 1%, above OSW+Percolator at 0.5%, 1% and 2%.
- **Calibration criterion: FAIL.** Bulk-null q ≤ 0.05 has FDP_c 0.068–0.072 (bar 0.06), and plain q is at 0.091. q ≤ 0.01 is calibrated (0.010).
- **Weights:** library dot product (−30) and spectral angle (+32) are nearly collinear, since the angle is a monotone function of the dot product. Under ridge 1e-3 they cancel. They are stable across folds, but redundant.
- **Cost:** 13:28 wall time against 12:58; peak RSS 25.1 GB.

Reading:
- The shape features are where the remaining raw 01 ID gap was. They make the score more discriminative.
- They also widen the target-versus-decoy null mismatch at 5%, which is the tail problem documented in CALIB. Bulk-null weighting corrects only part of it.
- Per the plan, this is not frozen. A post-hoc raw 02 run (`../../astral-raw02-shape-diagnostic-20261109`) checks whether the calibration failure replicates. It is a diagnostic, not a replication gate.

## Subgroup ablation (`shape_ablation.py`; exact tool-LDA replica, raw 01)

| arm (base = byion) | N at FDP_c 1% | plain q≤0.05 (FDP_c) | bulk-null q≤0.01 / 0.05 (FDP_c) |
|---|---|---|---|
| base | 1010 | 0.065 | 939 (0.006) / 1540 (0.050) |
| **+xcorr** | **1121 (+11%)** | 0.063 | 970 (0.008) / 1444 (**0.046**) |
| +library | 1222 | 0.093 | 1306 (**0.017**) / 1785 (**0.070**) |
| +signal | 1185 | 0.078 | 1123 (0.007) / 1636 (**0.065**) |
| +all (= tool) | 1233 | 0.091 | 1196 (0.010) / 1721 (**0.068**) |

- The library-agreement subgroup (dot product, Manhattan, RMSD, spectral angle) brings IDs and breaks calibration, even at q ≤ 0.01 (0.017).
- This is consistent with CALIB section 8. Entrapment winners are more library-like than decoys, so a feature that rewards library agreement separates entrapment (false) targets from decoys. That is exactly the false-target vs decoy null mismatch.
- The xcorr subgroup is the only one that adds IDs and keeps bulk-null q calibrated.
