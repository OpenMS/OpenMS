# RESULT_RTNULL: RT-refinement knockout does not fix calibration

26 November 2026. Frozen plan: PRECHECK_RTNULL.md, pushed as `83eed8ac1a` before jobs 44052409-44052418. All 10 cells exited 0. The fetched files match the remote SHA256 list and the baseline tables match `rtnull/baseline.sha256`. The unchanged frozen evaluator was used.

**Verdict: FAIL.** Neither set reaches a combined bound <= 1.00 %, and both lose IDs.

| Set | Combined G0 | Combined NORT | Median N@FDP_c 1 % G0 | Median N@FDP_c 1 % NORT | ID ratio |
|---|---:|---:|---:|---:|---:|
| D | 1.07 % (1.02-1.13) | 1.16 % (1.10-1.22) | 21,805 | 20,770 | 0.953x |
| H | 1.28 % (1.23-1.32) | 1.23 % (1.19-1.28) | 37,897 | 35,799 | 0.945x |

Ranges are the evaluator's 90 % Poisson intervals. Turning RT refinement off makes D worse and improves H by only 0.05 pp, at a cost of 5 % of IDs. P(entrapment beats its own decoy) is unchanged (D 0.536 to 0.543, H 0.533 to 0.530). The in-sample RT-centre asymmetry is therefore not a material source of the excess. Refinement stays on by default.

This knockout removes the whole second pass, not only the in-sample fit. A cross-fitted RT-centre variant is not worth building on this evidence. The remaining lead is the GBM fold-selection dependence in DIAG_CALIBRATION.md. Its FDP effect is not yet measured.

The first local fetch failed because the disk was full. Old scratch intermediates were deleted, then the outputs were fetched again into a new directory and checksum-verified. No remote output was changed. Full evaluator output: `rtnull/EVAL_RTNULL.txt`.
