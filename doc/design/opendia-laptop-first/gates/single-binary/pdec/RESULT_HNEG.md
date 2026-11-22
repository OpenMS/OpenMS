# RESULT_HNEG: C++ hard-negative GBM rescoring passes on r ≈ 1.56 (22 Nov 2026)

The gate is `PRECHECK_HNEG.md`, committed in docs before the code existed or any job ran.

- **Code:** OpenMS `d8ca86edeb` adds `-gbm_rescore` and `-gbm_hard_fraction` on top of the sibling default `e2e65938ce`. I reviewed it, and 27/27 tests pass, including the new `OpenDIAGBM_test`.
- **GBM:** written in-house, histogram-based, leaf-wise, and deterministic across thread counts. On a synthetic 3.15 M × 49 problem at 8 threads it fits in 81 s.
- **Container:** `opendia-pd10-gbmopt.sif`, sha256 `9b5df05f0972648ddd0c5ac7d9da2a9bf911c6b0ca810582935dd2fd5d2bada8`.
- **Runs:** `hneg-20261122/runs`. All 12 exited 0.
- **Evaluator:** `eval_pdec.py`, sha256 `bb03387d…`. Full output is in `EVAL_HNEG.txt`.

**Default unchanged.** With no flag, `precursors.tsv` is byte-identical to `sibcpp-20261121/runs/<R>_sib_s<seed>` (the sibling default) on all 4 runs.

## Numbers (r = 1.1225)

| Run | Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Gain | Wall |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| D | cur | 0 | 24,110 | 179 | 0.74 % | 1.40 % | 23,428 | | 7:47 |
| D | cur | 1 | 24,413 | 182 | 0.75 % | 1.41 % | 23,392 | | 5:23 |
| D | gbt | 0 | 25,995 | 220 | 0.85 % | 1.60 % | 24,714 | +5.5 % | 9:03 |
| D | gbt | 1 | 25,778 | 214 | 0.83 % | 1.57 % | 24,325 | +4.0 % | 10:33 |
| D | gbm | 0 | 25,169 | 204 | 0.81 % | 1.53 % | 24,103 | +2.9 % | 6:20 |
| D | gbm | 1 | 25,212 | 198 | 0.79 % | 1.48 % | 24,214 | +3.5 % | 7:18 |
| H | cur | 0 | 37,320 | 245 | 0.66 % | 1.24 % | 36,845 | | 8:59 |
| H | cur | 1 | 37,055 | 243 | 0.66 % | 1.24 % | 36,610 | | 10:25 |
| H | gbt | 0 | 39,834 | 298 | 0.75 % | 1.41 % | 38,903 | +5.6 % | 9:30 |
| H | gbt | 1 | 39,914 | 309 | 0.77 % | 1.46 % | 38,878 | +6.2 % | 8:02 |
| H | gbm | 0 | 39,129 | 283 | 0.72 % | 1.37 % | 38,420 | +4.3 % | 11:36 |
| H | gbm | 1 | 39,331 | 291 | 0.74 % | 1.40 % | 38,499 | +5.2 % | 8:53 |

**Pooled lower / combined bounds and pair fairness:**

| Run | cur | gbt | gbm |
|---|---|---|---|
| D | 0.74 / 1.41 %, P 0.499 | 0.84 / 1.59 %, P 0.542 | 0.80 / 1.51 %, P 0.538 |
| H | 0.66 / 1.24 %, P 0.507 | 0.76 / 1.44 %, P 0.518 | 0.73 / 1.38 %, P 0.515 |

## Rule (gbt)

| Criterion | Bar | Observed | |
|---|---|---|---|
| Mean N @ FDP_c 1 % gain | ≥ +3 % on D and H | D +4.7 %, H +5.9 % | pass |
| Every seed gains | yes | +5.5, +4.0, +5.6, +6.2 % | pass |
| Mean lower bound change | ≤ +0.2 pp | D +0.10 pp, H +0.10 pp | pass |
| Median wall | ≤ 1.5× cur | D 1.49×, H 0.90× | pass |

**PASS.** Freeze PLAN_GBT, the r ≈ 7 run on D seeds 0-4.

## Notes

- **Hard negatives help.** All-decoy-winner negatives (gbm) give D +3.2 % and H +4.7 %. Hard negatives (gbt) add about 1.5 pp more on D and 1.2 pp on H. The C++ port is smaller than the Python Bt prediction (r7 D +12.6 % / H +7.4 %, measured at r7, not r ≈ 1.56). The r7 plan will show whether the gap holds.
- **Calibration costs.**
  - **q looser:** the tool's own q is looser by about 0.1 pp on the lower bound, so slightly more entrapments pass at q ≤ 0.01. This is within the bar, and N at true FDP (FDP_c) still rises.
  - **Pair fairness:** on D it moves from 0.499 to 0.542. The plan bar is 0.56.
- **Wall time on D.** D wall is 1.49×, right at the bar. The GBM runs in both scoring calls, the RT-refine first pass and the final pass. At r ≈ 7, prediction over about 9 M windows will cost more, and the plan's cost criterion will decide.
