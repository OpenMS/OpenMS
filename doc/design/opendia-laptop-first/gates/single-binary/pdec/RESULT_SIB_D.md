# RESULT_SIB_D: PLAN_SIB passes on D at r ≈ 7 (21 Nov 2026)

The plan is `PLAN_SIB.md`, frozen and pushed in docs `1d5a75bf06` before any r ≈ 7 job ran.

- **Runs:** `entlib-sib-20261121/runs/D_{cur,sib}_s{0..4}`. All 10 exited 0.
- **Container:** `opendia-pd10-sibf.sif`, sha256 `c67d7d9f…`.
- **Evaluator:** `eval_pdec.py`, sha256 `bb03387d…`, r = 7.0539. Full output is in `EVAL_PLAN_SIB_D.txt`.

**Default unchanged.** CUR `precursors.tsv` is byte-identical to `entlib-rtr-20261121/runs/D_rtr_s<seed>` (the adopted default) on all 5 seeds.

## Numbers

| Seed | CUR N @ q ≤ 1 % | CUR N @ FDP_c 1 % | SIB N @ q ≤ 1 % | SIB N @ FDP_c 1 % | Gain |
|---|---:|---:|---:|---:|---:|
| 0 | 16,034 | 15,465 | 20,167 | 19,269 | +24.6 % |
| 1 | 15,877 | 15,369 | 20,281 | 19,865 | +29.3 % |
| 2 | 15,475 | 14,899 | 20,256 | 19,421 | +30.4 % |
| 3 | 16,091 | 15,464 | 20,381 | 19,704 | +27.4 % |
| 4 | 15,662 | 15,196 | 20,775 | 20,031 | +31.8 % |
| **median** | | **15,369** | | **19,704** | **+28.2 %** |

## Criteria

| # | Criterion | Bar | CUR | SIB | |
|---|---|---|---:|---:|---|
| 1 | Median N @ FDP_c 1 % | ≥ 1.03 × CUR | 15,369 | 19,704 (1.282×) | pass |
| 2 | Pooled lower bound | ≤ CUR + 0.0020 | 1.12 % | 1.10 % | pass |
| 2 | Pooled combined bound | ≤ CUR + 0.0030 | 1.27 % | 1.25 % | pass |
| 3 | Pair fairness P | ≤ 0.56 | 0.527 | 0.529 | pass |
| 4 | Median wall | ≤ 1.5 × CUR | 17:18 | 18:47 (1.09×) | pass |
| 4 | Peak RSS | ≤ 64 GB | 34.6 GB | 34.8 GB | pass |

**PASS on D.** As frozen, H is submitted now and reported either way.

**Where D stands.** On this run DIA-NN reports 26,006 at FDP_c 1 % (r7 reference). OpenDIA moves from 59 % to 76 % of that, and the 5 SIB seeds span 19,269-20,031.
