# RESULT_RTR_D: PLAN_RTR on D (21 Nov 2026)

## Verdict

**PASS, all four criteria, criterion 1 narrowly.** `-rt_refine` raises the median N at FDP_c 1 % from 14,895 to 15,369 (+3.2 %), against a frozen bar of +3.0 % (15,342). Unlike the earlier feature blocks, **the own-q calibration improves**: the pooled lower bound falls from 1.20 % to 1.12 %, and the combined bound from 1.37 % to 1.27 %.

- **Evaluation:** `pdec/eval_pdec.py` (sha256 `bb03387d…`), container `opendia-pd10-rtrefine2.sif` (sha256 `5e532efe…`) for both arms.
- **Default byte check:** the CUR arm's `precursors.tsv` is byte-identical to the PLAN_MT `mt` arm on all 5 seeds. That arm came from a different container (source 34c36861f2 with `-matched_filter`, which is today's default), so the flag-off path is unchanged at r ≈ 7.

## Per seed

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Wall (s) | Peak RSS (GB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CUR | 0 | 15,320 | 179 | 1.17 % | 1.33 % | 14,875 | 575 | 34.6 |
| CUR | 1 | 15,732 | 192 | 1.22 % | 1.39 % | 14,993 | 623 | 34.5 |
| CUR | 2 | 15,423 | 190 | 1.23 % | 1.41 % | 14,732 | 581 | 34.1 |
| CUR | 3 | 15,670 | 191 | 1.22 % | 1.39 % | 14,994 | 576 | 34.3 |
| CUR | 4 | 15,363 | 179 | 1.17 % | 1.33 % | 14,895 | 578 | 34.3 |
| RTR | 0 | 16,034 | 179 | 1.12 % | 1.27 % | 15,465 | 1,030 | 34.4 |
| RTR | 1 | 15,877 | 170 | 1.07 % | 1.22 % | 15,369 | 1,055 | 34.2 |
| RTR | 2 | 15,475 | 172 | 1.11 % | 1.27 % | 14,899 | 1,099 | 34.7 |
| RTR | 3 | 16,091 | 192 | 1.19 % | 1.36 % | 15,464 | 1,078 | 34.3 |
| RTR | 4 | 15,662 | 170 | 1.09 % | 1.24 % | 15,196 | 1,016 | 34.5 |

**RT refit, per RTR seed:**
- Training size was 15,320-15,732 IDs, with alpha 100 throughout.
- The p95 RT residual fell from 22.6-22.7 s to 17.2-17.3 s.
- The RT window narrowed from ±41.4 s to ±31.4-31.7 s, and candidates fell by 9 % (seed 0: 11.14 M to 10.16 M).

## Criteria

| # | Criterion | Value | Limit | Result |
|---|---|---|---|---|
| 1 | Median N @ FDP_c 1 %, RTR / CUR | 15,369 / 14,895 = **1.032** | ≥ 1.03 | pass (narrow) |
| 2a | Pooled lower bound, RTR − CUR | 1.12 − 1.20 = **−0.08 pp** | ≤ +0.20 pp | pass |
| 2b | Pooled combined bound, RTR − CUR | 1.27 − 1.37 = **−0.10 pp** | ≤ +0.30 pp | pass |
| 3 | Pair fairness, RTR | **0.527** (n 8,270, ±0.009) | ≤ 0.56 | pass |
| 4 | Median wall, RTR / CUR | 1,055 / 578 = **1.83** | ≤ 2.2 | pass |
| 4 | Peak RSS, RTR | **34.7 GB** | ≤ 64 GB | pass |

## Reading

- **The gain is small and seed-dependent.** Seeds 0, 1, 3 and 4 gain +2 to +4 %. Seed 2 is flat (14,899 against 14,732-14,994 for CUR). The median clears the bar by 27 IDs, so the sensitivity verdict is narrow and should not be read as a robust +3 %.
- **Calibration is the clearer effect.** On every seed the lower bound is at or below the best CUR seed, and the pooled 90 % intervals barely overlap (CUR 1.14-1.27 %, RTR 1.05-1.18 %). A narrower, better-centred RT window gives false targets fewer chances to find a matching peak, so fewer accepted IDs are entrapments.
- **The cost is about 1.8× wall time** for the second full pass, with no change in memory.
- **Against DIA-NN on D** (26,006 at FDP_c 1 %), this is 59 %, up from 57 %.

## Next

H (22474.d) is submitted after this verdict is committed and is reported either way. H has the larger RT error, so the precheck expected the bigger effect there (−38 % p95 against −27 % on D). Whether `-rt_refine` becomes the default is decided after H, given the 1.8× cost.
