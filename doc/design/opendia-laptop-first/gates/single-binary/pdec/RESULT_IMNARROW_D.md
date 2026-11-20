# RESULT_IMNARROW_D: narrower 1/K0 window on the r ≈ 7 entrapment, design run D: PASS (20 Nov 2026)

**Plan.** `PLAN_IMNARROW.md` was frozen at 07a4245472 (2026-11-20 08:07:05 UTC) and its sha256 checks OK. The codon
run directory `entlib-imnarrow-20261120` and its `run.sbatch` were created at 08:20:28 UTC, after the freeze. There was
no earlier r ≈ 7 run with `-im_tolerance_scale`.

**Setup.**
- Binary: OpenMS 0942ffbdf3, container pd10 (sha256 fca04be9…, checked at submission).
- Library: the PLAN_ODP predicted library.
- Arms: CUR is the default, NEW adds `-im_tolerance_scale 0.6`.
- 5 seeds each, 8 threads.
- Evaluator `eval_pdec.py` (sha256 bb03387d…). Output is in `EVAL_IMNARROW_D.txt`, and wall and RSS from `/usr/bin/time -v` are in `COST_IMNARROW_D.txt`.

## Per seed

| seed | CUR N at q ≤ 0.01 / E / N at FDP_c 1 % | NEW N at q ≤ 0.01 / E / N at FDP_c 1 % |
|---|---|---|
| 0 | 11,905 / 149 / 11,258 | 12,432 / 149 / 11,740 |
| 1 | 11,993 / 145 / 11,265 | 12,639 / 152 / 12,124 |
| 2 | 11,868 / 121 / 11,542 | 12,482 / 146 / 11,907 |
| 3 | 11,986 / 147 / 11,332 | 12,636 / 160 / 11,815 |
| 4 | 11,823 / 131 / 11,418 | 12,398 / 136 / 11,727 |

Every NEW seed is above every CUR seed at FDP_c 1 %: the lowest NEW seed has 11,727 and the highest CUR seed 11,542.

## Criteria

| # | Criterion | Value | Limit | Result |
|---|---|---|---|---|
| 1 | NEW median N at FDP_c 1 % ≥ 1.03 × CUR | 11,815 vs 11,332 (**+4.3 %**) | ≥ 11,671.96 | **pass** |
| 2a | NEW pooled lower bound ≤ CUR + 0.0020 | 0.011871 (743 / 62,587) | ≤ 0.013632 (CUR 0.011632) | pass |
| 2b | NEW pooled combined bound ≤ CUR + 0.0025 | 0.013554 | ≤ 0.015781 (CUR 0.013281) | pass |
| 3 | NEW pair win ≤ 0.56 | 0.501 (n 6,520) | ≤ 0.56 | pass |
| 4 | NEW median wall ≤ 1.4 × CUR, peak RSS ≤ 64 GB | 9:27.95 vs 9:40.48 (0.98×), 38.6 GiB | ≤ 1.4×, ≤ 64 GB | pass |

**Verdict: PASS.** Per the plan, H (22474.d) is submitted now and reported either way. Making 0.6 the default is a
separate commit, with ctests and a byte check.

## Reading

- The gain at r ≈ 7 (+4.3 %) is about half the r ≈ 1.56 seed-0 precheck (+8.3 %). The larger library is denser, so the effect is smaller.
- Calibration is unchanged within noise. The lower bound is 1.16 % → 1.19 %, with 90 % intervals 1.09-1.24 % and 1.12-1.26 %, which overlap. The r ≈ 1.56 improvement in calibration did not carry over.
- Pair win improved, 0.515 → 0.501.
- No extra cost at r ≈ 7: wall 0.98×, RSS 38.0 → 40.5 GB peak.
- The residual FDP excess of about 1.1-1.2× (open item b9) is unchanged by this change.
