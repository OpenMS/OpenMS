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

## Addendum: how DIA-NN rates the IDs we reject (`diann_view.py`, `PRECHECK_DIANNVIEW.txt`)

DIA-NN report columns for its non-entrapment q <= 0.01 IDs, grouped by our outcome (C++ seed 0). A = we accept the
right peak; R = right peak, below our cut; Wrong = our winner is > 6 s from DIA-NN's RT.

| | raw01 A | raw01 R | raw02 A | raw02 R |
|---|---|---|---|---|
| n | 1652 | 553 | 1678 | 493 |
| inside DIA-NN's own N@FDP_c 1 % set | 94 % | 73 % | 99 % | 96 % |
| DIA-NN Q.Value median | 0.0003 | 0.0036 | 0.0003 | 0.0016 |
| DIA-NN PEP median | 0.0003 | 0.017 | 0.0003 | 0.010 |
| abundance percentile median | 0.61 | 0.26 | 0.60 | 0.25 |
| Ms1.Profile.Corr median | 0.81 | 0.57 | 0.82 | 0.52 |

- R is DIA-NN's own low-confidence, low-abundance tail. On raw01, 27 % of R lie past DIA-NN's true-1 % cut, so
  part of the raw01 "gap" at DIA-NN's q <= 0.01 is DIA-NN's own excess FDP (0.0144), not missed evidence.
- Fragment mass accuracy, the other DIA-NN-style evidence not in our 24 inputs, was tested before (voxel v5,
  `gates/voxel-scoring-model/v4/RESULT_V5.md`) and failed. MS1 profile correlation is the column that separates A from R
  most, and the MS1 block (PLAN_MS1) is the closed test of that.
- A scan-level diagnostic (option b) would study weak peaks that DIA-NN itself rates at PEP 0.01-0.02. It is not
  expected to show a new evidence type, so it is not pursued.
