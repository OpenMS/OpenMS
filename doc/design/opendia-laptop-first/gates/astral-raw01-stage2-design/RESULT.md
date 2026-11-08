# Stage-2 design run on raw 01: result

Only target and decoy labels were used. Code commit `014cb0dbf1`. The rule in `PLAN.md` was written before the run.

## Rule check
- **Run:** completed in 12:17 wall time, with 22.8 GB peak.
- **Weights:** every fold gave finite weights with a positive `cnn_z` weight (9.70, 9.87, 9.73).
- **q ≤ 0.01 count:** 333, at least the 321 required.
- **Decision:** PASS. The `-stage2` defaults are frozen unchanged (3 iterations, 0.15/0.05, ±4 s, ±20 ppm). Nothing was tuned.

## Unchanged columns
The CNN logits, folds and RT-window flags are byte-identical to `../pipeline/crossfit_it0` (`cut -f1-6,9 | cmp`). Stage 2 changes only the winner choice and the q-values.

## Pooled paired q (experimental decoy-count estimate, not calibrated FDR)

| q ≤ | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 |
|---|---|---|---|---|---|
| CNN z only (it0) | 254 | 321 | n/a | 634 | n/a |
| stage 2 | 319 | 333 | 524 | 770 | 1013 |

- Per fold at q ≤ 0.01: 130 / 107 / 96, against 120 / 105 / 96 for CNN-only.
- Stage-2 ranking from the top: 200 targets with 0 decoys, 332 with 1, 396 with 4, 491 with 9, 676 with 24, 925 with 75.
- The q ≤ 0.01 count sits on a cliff. Because of the +1 in (D+1)/T, 4 decoys among 396 targets already gives q 0.0126. The gain is small at 1% and +21% at 5%.

## Per-fold LDA weights (standardised features)

| fold | cnn_z | frag_frac | frag_coel | lib_corr | ms2 ppm | ms1 ppm | ms1_coel | ms1_iso |
|---|---|---|---|---|---|---|---|---|
| 0 | 9.70 | −5.68 | 1.78 | 1.09 | −1.40 | −1.19 | 1.82 | −0.64 |
| 1 | 9.87 | −5.86 | 1.87 | 0.70 | −1.50 | −1.16 | 1.81 | −0.71 |
| 2 | 9.73 | −5.86 | 1.71 | 1.49 | −1.44 | −1.14 | 1.80 | −0.62 |

- **Stability:** the weights are stable across folds.
- **Negative weights:** `fragment_fraction` and `ms1_isotope` get negative weights. `fragment_fraction` correlates with `cnn_z` at r = 0.73, so the LDA uses it as a correction to `cnn_z`, not as evidence on its own. This is noted but not acted on, because changing it would mean tuning on this run.

## Feature medians

| group | n | frag_coel | lib_corr | ms1 ppm | ms1_coel |
|---|---|---|---|---|---|
| accepted targets | 333 | 0.94 | 0.95 | 0.46 | 0.84 |
| decoy winners with z > 4 | 7936 | 0.12 | 0.18 | 8.2 | 0.0 |
| target winners with z > 4 | 11071 | 0.17 | 0.33 | 7.0 | 0.13 |

## Limits
- Same-run design on one technical replicate.
- The feasibility rule had predicted about 690 targets above z 6. The fitted LDA gives fewer at the 1% cut, so the feasibility rule was optimistic.
