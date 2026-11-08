# MS1 isotope features: raw 01 result

Target/decoy labels only. `PLAN.md` was written before the run. Binary: commit `48e40f4d35` (`binary.sha256`).

- **Default unchanged:** `-stage2` with the isotope-collecting binary gives crossfit.tsv and features.tsv byte-identical to `../stage2-design` (`CMP_default.txt`).
- **Counts (pooled paired q, experimental decoy-count estimate):**

  | q ≤ | 0.005 | 0.01 | 0.02 | 0.05 |
  |---|---|---|---|---|
  | search | 331 | 717 | 1015 | 1412 |
  | search+isotope | 402 | 650 | 910 | 1354 |

- **Decision:** FAIL under the predeclared rule (at least 2% more at q ≤ 0.01). The isotope group stays implemented and off.
- **Weights:** stable and in the expected direction in all folds.
  - isotope correlation 0.79 to 0.91;
  - M−1 overlap −0.38 to −0.46;
  - isotope co-elution 0.70 to 0.75.
- **Medians by CNN z band** (target / decoy):

  | z band | isotope correlation | co-elution | M−1 overlap |
  |---|---|---|---|
  | z > 6 | 0.73 / 0.26 | 0.47 / 0.26 | 0.38 / 0.92 |
  | z 4 to 6 (where q ≤ 0.01 is decided) | 0.07 / 0.03 | 0.15 / 0.13 | 0.96 / 1.08 |

  - At z 4 to 6 the MS1 envelope is near noise for both classes. The features help at the very top (q ≤ 0.005: +71) and add noise at the 1% cut.
- **Not tuned:** wider core windows, more isotopes, or summed-intensity thresholds were not tried. Trying them on this run would be tuning.
