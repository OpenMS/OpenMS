# Stage-2 boundary group: raw 01 design diagnostic (written before the run)

- Raw 01 is design data and has already been used by the scoring ablation and four TOPAZ arms. This run is **descriptive and not a gate**. It only decides whether the boundary approach deserves a predeclared gate on an unused run.
- **One variable:** `-stage2_feature_set search` (frozen, 717 / 1015 / 1412) against `search+boundary`. Same binary, windows, seeds, CNN settings and seed.
- **Boundary group:** the window is re-picked with `OpenDIAPeakPicker` on its scan bins. Five features: log FWHM, log base width, co-elution over the picked peak, library correlation over the picked peak, and fragment apex agreement.
- **Reported:**
  - targets at pooled paired q ≤ 0.005 / 0.01 / 0.02 / 0.05;
  - per-fold weights;
  - medians of the boundary features for targets and decoys in CNN z bands z ≤ 0 and 4 < z ≤ 6. The bands are mostly null, so the features should not differ between targets and decoys there.
- **Reading, fixed now:**
  - a gain of at least 5% at q ≤ 0.01 (≥ 753) together with no target/decoy split in the null bands would justify a gate on an unused run;
  - anything less stops the rule-boundary approach. The learned segmenter is ported only if the gain is concentrated where the rule fails (low SNR).
- **Default equivalence:** `cmp` of `search/crossfit.tsv` and `features.tsv` against `../stage2-search-design/search`.
