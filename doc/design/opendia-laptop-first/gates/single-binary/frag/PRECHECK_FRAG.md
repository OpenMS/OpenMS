# PRECHECK_FRAG: run-specific fragment choice (diagnostic, no plan)

Question: would choosing which 6 of the 12 predicted fragments to search, based on what this run observes,
raise the evidence for DIA-NN-supported weak right-peak IDs (R) more than for score-matched decoys (D)?

Method (`frag_precheck.py`, output `PRECHECK_FRAG.txt`). Arm-E raw01/raw02 12-fragment arrays, C++ scorer seed 0.
- Observed rate: share of windows with central signal (bins 12-19 > 0) per library rank.
- Label-free, cross-fitted re-selection.
  - For each fold, ion-class observation rates (type x charge x ordinal bin) are learned from the other folds'
    target winners at q <= 0.01.
  - Per precursor, the 6 fragments with the highest rate x predicted intensity are kept.
- Compared: "support" (how many of the chosen 6 show central signal) for R vs D, A vs D, and E vs ED (specificity).

Results

| | raw01 orig | raw01 re-selected | raw02 orig | raw02 re-selected |
|---|---|---|---|---|
| R vs D support AUC | 0.708 | 0.710 | 0.725 | 0.721 |
| E vs ED support AUC | 0.502 | 0.489 | 0.500 | 0.514 |
| A vs D support AUC | 0.842 | 0.842 | 0.838 | 0.839 |

- Re-selection changes the set for 52 % of windows, but only about 1.15 fragments each.
- The searched 6 already show signal in 5.43 / 5.52 of 6 for R (5.90 / 5.91 for A).
- Observation falls smoothly with library rank for every group. No ion class is systematically missing in this run
  that re-selection could fix.
- The R deficit relative to A is spread over all ranks, which points to low abundance, not to the wrong fragments.
- Scoring already uses all 12 fragments (`library12_full`). Fragment choice only matters for the 6-fragment peak
  search, and PRECHECK_PEAK showed only 37-47 per raw lost there.

Verdict: no lever. Option (a) (run-specific fragment selection) is closed without a plan.
