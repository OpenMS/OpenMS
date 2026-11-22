# RESULT_GBT_D: PLAN_GBT passes on D at r ≈ 7 (22 Nov 2026)

The plan is `PLAN_GBT.md`, frozen and pushed in docs `9e54ff5e2f` before any r ≈ 7 job ran.

- **Runs:** `entlib-gbt-20261122/runs/D_{cur,gbt}_s{0..4}`. All 10 exited 0.
- **Container:** `opendia-pd10-gbmopt.sif`, sha256 `9b5df05f…`.
- **Evaluator:** `eval_pdec.py`, sha256 `bb03387d…`, r = 7.0539. Full output is in `EVAL_PLAN_GBT_D.txt`.

**Default unchanged.** CUR `precursors.tsv` is byte-identical to `entlib-sib-20261121/runs/D_sib_s<seed>` (the adopted sibling default) on all 5 seeds.

## Numbers

| Seed | CUR N @ q ≤ 1 % | CUR N @ FDP_c 1 % | GBT N @ q ≤ 1 % | GBT N @ FDP_c 1 % | Gain |
|---|---:|---:|---:|---:|---:|
| 0 | 20,167 | 19,269 | 22,204 | 21,468 | +11.4 % |
| 1 | 20,281 | 19,865 | 21,804 | 21,016 | +5.8 % |
| 2 | 20,256 | 19,421 | 22,153 | 21,396 | +10.2 % |
| 3 | 20,381 | 19,704 | 21,732 | 21,054 | +6.9 % |
| 4 | 20,775 | 20,031 | 22,017 | 21,459 | +7.1 % |
| **median** | | **19,704** | | **21,396** | **+8.6 %** |

## Criteria

| # | Criterion | Bar | CUR | GBT | |
|---|---|---|---:|---:|---|
| 1 | Median N @ FDP_c 1 % | ≥ 1.03 × CUR | 19,704 | 21,396 (1.086×) | pass |
| 2 | Pooled lower bound | ≤ CUR + 0.0020 | 1.10 % | 1.11 % | pass |
| 2 | Pooled combined bound | ≤ CUR + 0.0030 | 1.25 % | 1.27 % | pass |
| 3 | Pair fairness P | ≤ 0.56 | 0.529 | 0.544 | pass |
| 4 | Median wall | ≤ 1.5 × CUR | 19:06 | 19:55 (1.04×) | pass |
| 4 | Peak RSS | ≤ 64 GB | 34.6 GB | 34.8 GB | pass |

**PASS on D.** As frozen, H is submitted now and reported either way.

**Where D stands.** OpenDIA moves from 76 % to **82 % of DIA-NN** (26,006), and every GBT seed is above every CUR seed. At r ≈ 7 the gain (+8.6 %) is larger than at r ≈ 1.56 (+4.7 %), but smaller than the Python Bt estimate (+12.6 %). The r ≈ 7 GBM cost is about 50 s per run, small next to extraction.
