# RESULT_EXTRA_H: PLAN_EXTRA held-out run H (19 Nov 2026, reported only)

**Order.** The D verdict was committed in d44389e9bd at 18:38:56Z, and the H jobs (36504048-36504069) were submitted at
18:39:06Z. Setup and evaluator are as on D. Output is in `EVAL_EXTRA_H.txt`. CUR matches the earlier PLAN_ODP H CUR
numbers (RESULT_H: about 20.7 k IDs at FDP_c 1 %).

| | CUR (default) | NEW (`-extra_scores`) |
|---|---|---|
| N at own q ≤ 0.01, per seed | 21,386 / 21,415 / 21,232 / 21,342 / 21,309 | 26,224 / 26,317 / 26,357 / 26,311 / 26,246 |
| **median N at FDP_c 1 %** (range) | **20,745** (20,617-20,810) | **25,643** (25,578-25,759) |
| pooled lower bound at q ≤ 0.01 (90 %) | 0.0106 (0.0101-0.0111) | 0.0107 (0.0103-0.0112) |
| pooled combined bound at q ≤ 0.01 (90 %) | 0.0121 (0.0115-0.0127) | 0.0123 (0.0117-0.0128) |
| P(entrapment beats own decoy, best q ≤ 0.05) | 0.505 | 0.483 |
| wall / peak RSS per run | 10:34-11:06 / 32.5-33.3 GB | 10:46-13:32 / 32.7-34.4 GB |

**Applying the D criteria to H (for information; the plan's verdict is on D):**
- N at FDP_c 1 % is 1.236× CUR, against a limit of 1.10.
- The lower bound rises by 0.0001 and the combined bound by 0.0002, against limits of 0.0025 and 0.0030.
- Pair win is 0.483, against a limit of 0.56.

All of these hold. Every NEW seed is above every CUR seed.

## Summary of PLAN_EXTRA
| | D | H |
|---|---|---|
| IDs at 1 % true FDP, CUR to NEW | 8,233 → 11,332 (+37.6 %) | 20,745 → 25,643 (+23.6 %) |
| share of DIA-NN 2.2 (26,006 / 41,272) | 32 % → 44 % | 50 % → 62 % |
| own q ≤ 0.01 true FDP, lower / combined | 1.06 / 1.21 % → 1.16 / 1.33 % | 1.06 / 1.21 % → 1.07 / 1.23 % |

- On H the gain comes at no calibration cost.
- On D calibration loosens a little but stays within the frozen tolerance.
- Absolute calibration is still the open PLAN_ODP item for both arms: the lower bound is about 1.06 % against the
  1.0 % limit.

**Next, per the plan's outcome:** propose `-extra_scores` as the default in a separate commit, with ctests and a byte
check that the new default equals the tested NEW arm.
