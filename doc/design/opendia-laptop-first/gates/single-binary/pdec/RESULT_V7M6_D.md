# RESULT_V7M6_D: PLAN_V7M6 fails on D at r ≈ 7, on wall time only (23 Nov 2026)

The plan is `PLAN_V7M6.md`, frozen and pushed in docs `b7c772816c` before any plan job ran (sha256 `ad3dc855`).

- **Runs:** `entlib-v7m6-20261123/runs/D_{cur,v7m6}_s{0..4}`. All 10 exited 0.
- **Container:** `opendia-pd10-fastint.sif`, sha256 `8de3a968…`.
- **Evaluator:** `eval_pdec.py`, sha256 `bb03387d…`. Full output is in `EVAL_PLAN_V7M6_D.txt`.

**Default unchanged.** CUR `precursors.tsv` is byte-identical to `entlib-gbt-20261122/runs/D_gbt_s<seed>` (the adopted GBM default) on all 5 seeds.

## Numbers

| Seed | CUR N @ FDP_c 1 % | V7M6 N @ FDP_c 1 % | Gain | CUR wall | V7M6 wall |
|---|---:|---:|---:|---:|---:|
| 0 | 21,468 | 24,324 | +13.3 % | 20:16 | 26:39 |
| 1 | 21,016 | 24,249 | +15.4 % | 18:55 | 26:50 |
| 2 | 21,396 | 23,896 | +11.7 % | 18:57 | 30:09 |
| 3 | 21,054 | 24,551 | +16.6 % | 18:39 | 28:43 |
| 4 | 21,459 | 23,962 | +11.7 % | 20:24 | 35:28 |
| **median** | **21,396** | **24,249** | **+13.3 %** | **18:57** | **28:43** |

## Criteria

| # | Criterion | Bar | CUR | V7M6 | |
|---|---|---|---:|---:|---|
| 1 | Median N @ FDP_c 1 % | ≥ 1.03 × CUR | 21,396 | 24,249 (1.133×) | pass |
| 2 | Pooled lower bound | ≤ CUR + 0.0020 | 1.11 % | 1.11 % | pass |
| 2 | Pooled combined bound | ≤ CUR + 0.0030 | 1.27 % | 1.27 % | pass |
| 3 | Pair fairness P | ≤ 0.56 | 0.544 | 0.540 | pass |
| 4 | Median wall | ≤ 1.5 × CUR | 18:57 | 28:43 (**1.515×**) | **fail** |
| 4 | Peak RSS | ≤ 64 GB | 34.4 GB | 47.1 GB | pass |

**Verdict: FAIL on criterion 4 (median wall 1.515× against 1.5×), with criterion 1 passing.** Per the plan, the gain is reported with its cost and the defaults stay. Nothing is retuned, and H is not run as a plan step.

## Reading

- **The gain is real and robust.** Every seed gains 11.7-16.6 %, and the q calibration and pair fairness are unchanged or slightly better.
- **The wall bar is missed by 1 %, on unmatched nodes.** The plan runs on standard codon nodes, as PLAN_GBT did. The ten jobs landed on six AMD nodes. The slowest V7M6 run (seed 4, 35:28) was the only one on `hl-codon-cn-002`, while the other four V7M6 runs took 26:39-30:09. On matched Intel nodes (RESULT_V7M6T2) the same comparison was 1.375×. By the frozen rule this is still a fail. I do not re-time, drop seeds or move the bar.
- **What V7M6 would bring:** D at r ≈ 7 from 21,396 to 24,249 at true FDP 1 %, which is 93 % of DIA-NN's 26,006.

## Options (need a decision)

1. **Accept the fail.** V7M6 stays an opt-in setting (`-valley_fraction 0.7 -max_peaks 6`).
2. **Further output-identical speed-ups,** then a new frozen plan. After the resolver change, the remaining cost is extraction and scoring of 1.39× more windows.
3. **Change the timing rule for future plans** (for example, matched nodes or ABBA within a job). That is a rule change and needs your call. It would not rescue this plan.
