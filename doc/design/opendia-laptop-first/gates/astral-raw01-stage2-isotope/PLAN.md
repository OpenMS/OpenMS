# MS1 isotope features: raw 01 design run (written before the run)

Target/decoy labels only.
- **Code change:** `TraceCollector` evidence now also records MS1 M−1/z and M+0..M+3 intensities in the core window, plus M+1 and M+2 traces. They feed three opt-in features in the stage-2 `isotope` group:
  - averagine correlation of M+0..M+3, the same model OpenSWATH's MS1 isotope score uses (OpenSwathWorkflow defaults to `-ms1_isotopes 3`);
  - M−1 overlap ratio;
  - M+1/M+2 co-elution with M+0.
- **Runs, same binary and same frozen CNN flags:**
  1. `-stage2` default (`extraction`): must be byte-identical to `../stage2-design` (crossfit.tsv and features.tsv). This checks that the isotope collection changed nothing.
  2. `-stage2_feature_set search+isotope -stage2_search_features ../pipeline/search.features.tsv`.
- **Comparison:** run 2 against `../stage2-search-design/search`, which has the same features minus the isotope group.
- **Decision rule, fixed now:**
  - The isotope group is adopted if run 2 gives at least 2% more targets at q ≤ 0.01 than the `search` run, and every fold has finite weights with positive `cnn_z`.
  - Otherwise it stays implemented but off, and the result is recorded.
  - No parameter is tuned on this run.
