# Native stage-2 with search peak features: raw 01 result

Target/decoy labels only. `PLAN.md` was written before the run.

- **Default unchanged:** `-stage2` with the refactored binary gives crossfit.tsv and features.tsv byte-identical to `../stage2-design` (`CMP_default.txt`).
- **`-stage2_feature_set search`:** CNN logit, transition support, log matched scans, log peak width and |library RT error|.
  - 717 / 1015 / 1412 targets at pooled paired q ≤ 0.01 / 0.02 / 0.05. The offline replica gave 717 / 1014 / 1412, so it reproduces.
  - Per fold at q ≤ 0.01: 254 / 238 / 225.
  - The CNN logit, fold and window columns are byte-identical (`CMP_search.txt`).
- **Weights (standardised):** stable across folds, with `cnn_z` positive in every fold.

  | | cnn_z | support | log scans | log width | abs ΔRT |
  |---|---|---|---|---|---|
  | range over folds | 5.22–5.35 | 6.19–6.32 | −7.56 to −8.01 | 3.33–3.70 | −2.30 to −2.40 |

  - Matched scans gets a negative weight once width is in the model. A target needs fewer scans at a given width, meaning a sharper peak.
- **Decision:** PASS. The search set is frozen as the stage-2 design candidate.
- **Scale:** 717 against 333 for the extraction set and 321 for CNN-only. For reference, OpenSWATH+Percolator on the same library gives 1441 on raw 01.
- **Limits:**
  - These are experimental decoy-count estimates, not entrapment-checked.
  - It is the same run the design was made on.
- **Resources:** 14:19 wall time, 23.2 GB peak.
