# Why we still lag OSW+Percolator on raw 01 (diagnostic, design run only)

This is a post-hoc diagnostic on raw 01, which is a heavily reused design run. It is not a gate. The q-values are paired-pooled decoy-count estimates computed offline in Python (`lag.py`). They are not the tool's output.

## Funnel (search+boundary, q <= 0.01)

Of the 1462 OSW 1% IDs on raw 01:

- 1379 have our precursor winner at the OSW RT;
- 710 of those are accepted;
- 672 are right peaks that are not accepted. Their median q is 0.032.

So peak finding is not the limit. Ranking right peaks above decoys is.

## Is the information there? (per-feature AUC)

This compares the 604 missed right-peak targets against the 1146 decoy winners in the same stage-2 score band.

- **Strong separation remains:**
  - `transition_support` 0.84
  - `library_correlation` 0.82
  - `boundary_coelution` and `boundary_library_correlation` 0.81
  - `boundary_apex_agreement` 0.80
  - `ms1_abs_ppm` 0.26 and `ms2_abs_ppm` 0.32 (that is, 0.74 and 0.68 when inverted)
  - `abs_library_rt_error` 0.34
- **Weak separation:** `cnn_z` 0.57.

The current linear stage-2 combination leaves strong single-feature signal unused inside the band.

## Combiner test (cross-fitted, semi-supervised, same folds)

| Combiner / features | q <= 0.005 | 0.01 | 0.02 | 0.05 | OSW IDs in 1% | R |
|---|---|---|---|---|---|---|
| LDA search+boundary (≈ tool) | 361 | 745 | 1001 | 1440 | 709 | 0.968 |
| MLP 32-16 search+boundary, seed 0 | 283 | 947 | 1177 | 1536 | 849 | 0.973 |
| MLP seeds 1 / 2 / 3 (0.01 only) | | 1131 / 761 / 972 | | 1536 / 1531 / 1550 | 987 / 667 / 875 | |
| MLP 5-seed ensemble A | 309 | 1009 | 1258 | 1534 | 919 | 0.984 |
| MLP 5-seed ensemble B | 353 | 935 | 1174 | 1544 | 858 | 0.980 |
| LDA +extraction | 473 | 641 | 1111 | 1351 | 619 | 0.968 |
| MLP +extraction, seed 0 | 0 | 663 | 1039 | 1441 | 624 | 0.955 |

## Reading

1. **The combiner is a real limit.** A small non-linear combiner on the same 10 features adds about 7% at q <= 0.05, and this holds across all seeds (1531–1550 against 1440). At q <= 0.01 it adds about 25–35% (ensembles 935–1009 against 745). The OSW 1% IDs it recovers rise from 709 to 858–919.
2. **The top of the list is unstable.** For single seeds the 1% count ranges from 761 to 1131, and the 0.5% count collapses on some seeds. The q <= 0.01 cut sits in a thin region where a few decoys move the count a lot. Ensembling narrows the range but does not remove it. A native combiner would therefore need a seed ensemble or a smoother fit, plus a predeclared check of variance across seeds.
3. **Adding the extraction features hurts both combiners.** The cause was not investigated. One hypothesis is that the semi-supervised loop locks onto the wrong structure early (R drops for the MLP). The extraction features are informative in isolation (AUC 0.77–0.82 for coelution and library correlation), so the loss is not explained by a lack of information.
4. **Null symmetry looks acceptable.** R at the decoy median is 0.97–0.98 for all rows. That supports the counts being honest decoy-count estimates rather than overfit. It is not entrapment-validated.
5. **CNN discrimination is weak in the miss band** (AUC 0.57). The CNN mostly decides which peak is chosen, not whether that peak is a target.

## Caveats

- This is a single run, and raw 01 is reused, so the headline numbers are optimistic.
- The Python q computation is paired-pooled and only approximates the tool's.
- OSW overlap is a diagnostic only. It is never used as a label.

## Implication

The next lever is a non-linear stage-2 combiner (a small MLP ensemble, or gradient-boosted stumps as PyProphet XGBoost does), opt-in, with defaults unchanged. It must be evaluated with a predeclared rule on an unused design run (not raw 01, 02 or 03) and include seed variance and an entrapment check.
