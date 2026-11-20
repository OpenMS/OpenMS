# RESULT_MT_D: PLAN_MT on D (20 Nov 2026)

## Verdict

**PASS, all four criteria.** On the r ≈ 7 entrapment, `-matched_filter` raises the median N at FDP_c 1 % from 13,985 to 14,895 (+6.5 %), and every MT seed is above every CUR seed. The own-q bounds rise by about 0.1 percentage points, inside the frozen limits. Cost is unchanged.

- **Evaluation:** `pdec/eval_pdec.py` (sha256 `bb03387d…`). The container was `opendia-pd10-34c36861f2.sif` (sha256 `1cad85f5…`), for both arms.
- **Default byte check:** the CUR arm's `precursors.tsv` is byte-identical to the PLAN_FAINT CUR arm on all 5 seeds, which comes from a different container built from source `2257c55c7a`. So the flag-off path is unchanged at r ≈ 7 as well.

## Per seed

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Wall (s) | Peak RSS (GB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CUR | 0 | 14,562 | 157 | 1.08 % | 1.23 % | 13,924 | 596 | 34.0 |
| CUR | 1 | 14,670 | 162 | 1.10 % | 1.26 % | 14,267 | 578 | 34.1 |
| CUR | 2 | 14,612 | 171 | 1.17 % | 1.34 % | 14,009 | 561 | 34.0 |
| CUR | 3 | 14,578 | 159 | 1.09 % | 1.25 % | 13,985 | 715 | 34.1 |
| CUR | 4 | 14,390 | 155 | 1.08 % | 1.23 % | 13,889 | 772 | 34.1 |
| MT | 0 | 15,320 | 179 | 1.17 % | 1.33 % | 14,875 | 568 | 34.2 |
| MT | 1 | 15,732 | 192 | 1.22 % | 1.39 % | 14,993 | 636 | 34.4 |
| MT | 2 | 15,423 | 190 | 1.23 % | 1.41 % | 14,732 | 589 | 34.4 |
| MT | 3 | 15,670 | 191 | 1.22 % | 1.39 % | 14,994 | 578 | 34.3 |
| MT | 4 | 15,363 | 179 | 1.17 % | 1.33 % | 14,895 | 551 | 34.4 |

## Criteria

| # | Criterion | Value | Limit | Result |
|---|---|---|---|---|
| 1 | Median N @ FDP_c 1 %, MT / CUR | 14,895 / 13,985 = **1.065** | ≥ 1.05 | pass |
| 2a | Pooled lower bound, MT − CUR | 1.201 − 1.104 = **+0.097 pp** | ≤ +0.20 pp | pass |
| 2b | Pooled combined bound, MT − CUR | 1.371 − 1.261 = **+0.111 pp** | ≤ +0.30 pp | pass |
| 3 | Pair fairness, MT | **0.533** (n 8,262, ±0.009) | ≤ 0.56 | pass |
| 4 | Median wall, MT / CUR | 578 / 596 = **0.97** | ≤ 1.3 | pass |
| 4 | Peak RSS, MT | **34.4 GB** | ≤ 64 GB | pass |

## Reading

- **The gain shrank from r ≈ 1.56 to r ≈ 7, as it did for PLAN_FAINT, but it stays above the bar.** It went from +9.0 % to +6.5 %. PLAN_FAINT went from +6.3 % to +4.0 % and failed.
- **The calibration cost is real but small.**
  - The pooled lower bound goes from 1.10 % to 1.20 %, and 90 % intervals overlap (1.04-1.17 % and 1.14-1.27 %).
  - The tool accepts +4,696 at q ≤ 1 % (pooled, 5 seeds) with +127 entrapments. With k = 1.142 that is about 145 estimated false, so about 3.1 % of the additions are false. That is higher than the 1.3 % in the base set, so these additions are a weaker tail. At equal true FDP the tool still gains, which is what criterion 1 measures.
- **Pair fairness moved from 0.499 to 0.533,** which is inside the limit but the largest move of any recent block. A false target now beats its own decoy a little more often, consistent with the bound increase. It is watched on H.
- **Position against DIA-NN on D (26,006 at FDP_c 1 %):** this is 57 % (was 54 %).

## Next

H (22474.d) is submitted after this verdict is committed and is reported either way. Making `-matched_filter` the default is a separate commit with a byte check against the tested arm.
