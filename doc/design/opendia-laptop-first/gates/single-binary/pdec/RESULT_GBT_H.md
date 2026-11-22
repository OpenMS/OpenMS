# RESULT_GBT_H: PLAN_GBT on held-out H (22 Nov 2026, reported, no criteria)

H (22474.d) was submitted after the D verdict was committed in docs `1362c68d5a`. It used the same container, library, evaluator and arms as D. The CUR arm's `precursors.tsv` is byte-identical to `entlib-sib-20261121/runs/H_sib_s<seed>`, the adopted default, on all 5 seeds. All 10 runs exited 0. Full output is in `EVAL_PLAN_GBT_H.txt`.

## Per seed

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Gain | Wall |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CUR | 0 | 35,717 | 414 | 1.16 % | 1.32 % | 34,953 | | 22:19 |
| CUR | 1 | 35,445 | 386 | 1.09 % | 1.24 % | 34,839 | | 29:28 |
| CUR | 2 | 35,457 | 400 | 1.13 % | 1.29 % | 34,684 | | 31:13 |
| CUR | 3 | 35,595 | 405 | 1.14 % | 1.30 % | 34,989 | | 22:44 |
| CUR | 4 | 35,675 | 434 | 1.22 % | 1.39 % | 34,871 | | 27:28 |
| GBT | 0 | 38,507 | 498 | 1.29 % | 1.48 % | 37,339 | +6.8 % | 32:15 |
| GBT | 1 | 38,359 | 470 | 1.23 % | 1.40 % | 37,342 | +7.2 % | 25:25 |
| GBT | 2 | 38,462 | 501 | 1.30 % | 1.49 % | 37,349 | +7.7 % | 30:26 |
| GBT | 3 | 38,531 | 503 | 1.31 % | 1.49 % | 37,428 | +7.0 % | 30:57 |
| GBT | 4 | 38,380 | 504 | 1.31 % | 1.50 % | 37,213 | +6.7 % | 24:13 |

## Pooled

| Arm | Lower (90 %) | Combined (90 %) | Median N @ FDP_c 1 % | Pair win | Median wall | Peak RSS |
|---|---|---|---:|---:|---:|---:|
| CUR | 1.15 % (1.10-1.19) | 1.31 % (1.26-1.36) | 34,871 | 0.521 | 27:28 | 29.4 GB |
| GBT | 1.29 % (1.25-1.33) | 1.47 % (1.42-1.52) | 37,342 | 0.540 | 30:26 (1.11×) | 28.6 GB |

## Reading

- **Gain.** Median N at FDP_c 1 % rises by **+7.1 %**, and every GBT seed (37,213-37,428) is above every CUR seed (34,684-34,989). OpenDIA moves from 84 % to **90 % of DIA-NN** (41,272).
- **Calibration cost on H.**
  - The tool's own q loosens: the pooled lower bound rises by +0.14 pp and the combined by +0.16 pp. If PLAN_GBT's criterion 2 bars were applied to H, it would pass them (+0.20 / +0.30 pp).
  - The lower bound at own q ≤ 1 % is now 1.29 %, against 1.15 %. The extra IDs are real at true FDP 1 % (FDP_c), but the reported q is slightly optimistic on this run.
  - Pair fairness goes from 0.521 to 0.540, under the 0.56 bar.
- **Cost.** Wall time is 1.11× on node-noisy runs, and memory is unchanged.
- **Verdict.** D passed every PLAN_GBT criterion (RESULT_GBT_D) and H confirms the gain. Next is the default flip, in a separate OpenMS commit with a byte check both ways. The looser q on H is recorded as a known cost to watch: it is the first lever since the matched filter that moved the tool's own calibration in the optimistic direction.
