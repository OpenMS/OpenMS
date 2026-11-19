# RESULT_IM2D_H: PLAN_IM2D held-out run H (19 Nov 2026, reported only)

**Order.** The D verdict (FAIL, criterion 1) was committed in 95383ade2c at 21:23:10Z. The H jobs (36803215-36803233)
were submitted at 21:23:21Z. Setup and evaluator are as on D. Output is in `EVAL_IM2D_H.txt`. CUR seed 0
`precursors.tsv` is byte-identical to PLAN_EXTRA H NEW seed 0.

| | CUR (default) | NEW (`-im_peaks 8`) |
|---|---|---|
| N at own q ≤ 0.01, per seed | 26,224 / 26,317 / 26,357 / 26,311 / 26,246 | 26,872 / 27,076 / 26,778 / 27,111 / 27,067 |
| **median N at FDP_c 1 %** (range) | **25,643** (25,578-25,759) | **26,302** (26,187-26,422) |
| pooled lower bound at q ≤ 0.01 (90 %) | 0.0107 (0.0103-0.0112) | 0.0108 (0.0103-0.0113) |
| pooled combined bound at q ≤ 0.01 (90 %) | 0.0123 (0.0117-0.0128) | 0.0123 (0.0118-0.0129) |
| pair win | 0.483 | 0.489 |
| median wall / max peak RSS | 11:19 / 33.2 GB | 16:37 / 42.8 GB |

**D criteria applied to H, for information only (the verdict is on D and stands):**
- N at FDP_c 1 % is 1.026× CUR (limit 1.03), so it would fail narrowly.
- The lower bound rises by 0.0001 and the combined bound by 0.0000, so calibration passes.
- Pair win is 0.489, which passes.
- Wall is 1.47× CUR (limit 1.5) and RSS is 42.8 GB, which passes.

On H every NEW seed is above every CUR seed (+2.6 %). On D the gain is +0.8 %, with overlapping seed ranges.

## Summary of PLAN_IM2D
| | D | H |
|---|---|---|
| IDs at 1 % true FDP, CUR to NEW | 11,332 → 11,420 (+0.8 %) | 25,643 → 26,302 (+2.6 %) |
| own q ≤ 0.01 true FDP, lower / combined | 1.16 / 1.33 % → 1.19 / 1.36 % | 1.07 / 1.23 % → 1.08 / 1.23 % |
| wall | ×1.49 | ×1.47 |

- **RT × mobility picking is safe but small.** It does not hurt calibration or fairness, and it gives +1-3 % IDs at about
  +50 % run time.
- **It stays opt-in (`-im_peaks`).**
- **The cost is the extra candidates, not the picker.** A version that adds band peaks only where they can win (scored
  on the band itself) would be needed to change the trade-off. That needs its own precheck.
