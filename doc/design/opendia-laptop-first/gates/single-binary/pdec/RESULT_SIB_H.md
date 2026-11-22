# RESULT_SIB_H: PLAN_SIB on held-out H (22 Nov 2026, reported, no criteria)

H (22474.d) was submitted after the D verdict was committed in docs `896fbeff9f`. It used the same container, library, evaluator and arms as D. The CUR arm's `precursors.tsv` is byte-identical to `entlib-rtr-20261121/runs/H_rtr_s<seed>`, the adopted default, on all 5 seeds. All 10 runs exited 0. Full output is in `EVAL_PLAN_SIB_H.txt`.

## Per seed

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Gain | Wall | Peak RSS (GB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CUR | 0 | 31,260 | 375 | 1.20 % | 1.37 % | 30,543 | | 21:18 | 26.8 |
| CUR | 1 | 31,145 | 352 | 1.13 % | 1.29 % | 30,618 | | 22:26 | 27.2 |
| CUR | 2 | 31,387 | 374 | 1.19 % | 1.36 % | 30,665 | | 28:03 | 27.0 |
| CUR | 3 | 31,490 | 367 | 1.17 % | 1.33 % | 30,868 | | 20:25 | 26.6 |
| CUR | 4 | 31,282 | 359 | 1.15 % | 1.31 % | 30,552 | | 28:49 | 33.1 |
| SIB | 0 | 35,717 | 414 | 1.16 % | 1.32 % | 34,953 | +14.4 % | 22:32 | 28.3 |
| SIB | 1 | 35,445 | 386 | 1.09 % | 1.24 % | 34,839 | +13.8 % | 32:28 | 28.1 |
| SIB | 2 | 35,457 | 400 | 1.13 % | 1.29 % | 34,684 | +13.1 % | 25:07 | 28.9 |
| SIB | 3 | 35,595 | 405 | 1.14 % | 1.30 % | 34,989 | +13.4 % | 26:22 | 27.9 |
| SIB | 4 | 35,675 | 434 | 1.22 % | 1.39 % | 34,871 | +14.1 % | 25:24 | 27.8 |

## Pooled

| Arm | Lower (90 %) | Combined (90 %) | Median N @ FDP_c 1 % | Pair win | Median wall |
|---|---|---|---:|---:|---:|
| CUR | 1.17 % (1.12-1.21) | 1.33 % (1.28-1.38) | 30,618 | 0.517 | 22:26 |
| SIB | 1.15 % (1.10-1.19) | 1.31 % (1.26-1.36) | 34,871 | 0.521 | 25:24 (1.13×) |

## Reading

- **Gain.** Median N at FDP_c 1 % rises by +13.9 %, and every seed gains between +13.1 % and +14.4 %. This matches the r ≈ 1.56 precheck (+13.5 %), so the gain does not shrink on the larger library here.
- **Calibration and fairness.** The pooled bounds are slightly lower and pair fairness is unchanged.
- **Cost.** Wall time is up 1.13× and memory by about 1 GB.
- **Against DIA-NN.** DIA-NN reports 41,272 on this run at r ≈ 7, so OpenDIA moves from 74 % to 84 % of it.
- **Verdict.** D passed every PLAN_SIB criterion (RESULT_SIB_D) and H confirms it. Next is the default flip, in a separate OpenMS commit with a byte check both ways.

## Default flip (22 Nov 2026)

- **Commit:** OpenMS `e2e65938ce` makes the sibling rescoring the default and adds `-no_sibling_features`. 26/26 tests pass.
- **Container:** `opendia-pd10-sibdef.sif`, sha256 `78cc3a625e28a15c2bb5c30f2a42606ae8aac88dcd6ddf6b67f4a936e9474f7b`. Runs are in `sibdef-20261122`.
- **Byte check on r ≈ 1.56 D seed 0:**
  - The default run is byte-identical to the tested `-sibling_features` arm (`sibcpp-20261121/runs/D_sib_s0`) on all four outputs: `precursors.tsv`, the run table, `calibration.json` and `features.tsv`.
  - `-no_sibling_features` is byte-identical to the previous default (`D_cur_s0`) on all three outputs.
