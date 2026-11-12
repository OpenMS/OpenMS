# V13: precursor-aware full-spectrum deconvolution diagnostic

**Status: diagnostic only; no training.** Raw01, 12 November 2026.

## Question and frozen rule
Could a label-blind nonnegative fit of the complete MS2 spectrum separate D11FWT stage-D candidates from score- and XIC-area-matched decoy winners, while leaving entrapment targets near the decoy distribution? The frozen V13b rule required candidate-ablation AUC ≥0.60 for D vs DEC, ≤0.53 for ENT vs DEC, and at least 80% valid fits per group. This exploratory diagnostic does not change the D11FWT comparison base or its accepted-ID pass rules.

The frozen panel contains all 866 stage-D winners, with one score/log-area matched control per D from each of DEC and ENT. Controls were matched with replacement. The label-free feature file contains 1,871 unique precursor IDs. The extractor saw only precursor IDs/apex RTs, mzML spectra, calibration, and target+decoy library precursor/transition records. Group labels were joined only after extraction.

## V13b result

| Measure | D | DEC | ENT | D vs DEC AUC | ENT vs DEC AUC |
|---|---:|---:|---:|---:|---:|
| Candidate ablation SSE increase, median | 0.0006 | 0.0001 | 0.0002 | **0.6872** | **0.5348** |
| Candidate coefficient fraction, median | 0.0192 | 0.0082 | 0.0083 | 0.6751 | 0.5220 |
| Candidate matched centroids, median | 4 | 4 | 4 | 0.6187 | 0.4632 |
| Normalized full-fit residual, median | 0.9544 | 0.8972 | 0.9286 | 0.6211 | 0.5524 |

All 866 rows in each group had a valid fit. Median local library neighborhoods were large: 555 precursors for D, 640.5 DEC, and 561.5 ENT; median fitted template counts were 75, 120, and 102. The normalized residual remains high, so the library templates do not explain most of the full spectrum. The raw residual AUCs (0.621 D-vs-DEC, 0.552 ENT-vs-DEC) mean D had a larger residual, not a better fit. Paired on each D, D's ablation score exceeded its matched DEC in 71.1% of pairs and its matched ENT in 66.2%.

## Decision
**The frozen gate fails narrowly.** D-vs-DEC ablation AUC clears the minimum by 0.087, and the coefficient-fraction feature places ENT-vs-DEC near chance at 0.522. However, the primary candidate-ablation ENT-vs-DEC AUC is 0.5348, slightly above the predeclared maximum of 0.53. The diagnostic therefore does not authorize adding the feature or starting model training. This is a useful signal for D against decoys, but not yet sufficiently decoy/entrapment-symmetric under the frozen rule.

The very large and correlated local template sets and high residuals limit physical interpretation of NNLS coefficients. Treat the reported features as ranking diagnostics, not quantitative fragment ownership.

## Invalid V13 attempt and correction
The first extraction multiplied scan-index RT values by 60. Direct comparison of the same spectrum index in `scan_index.parquet` and `OnDiscMSExperiment.getRT()` showed both values are already in seconds. The V13 features/evaluation are invalid and excluded. V13b uses scan-index RT directly, adds a regression test for seconds-based scan selection, and writes a separate output. The V13b plan and script/panel hashes were frozen at `2026-11-12T00:03:05Z` before corrected extraction.

## Reproducibility

- Plan: `PLAN_V13b.md`, SHA-256 `f48d62d8182ff5c39f523db8d89f6452f5e974fc75a446550b2c12818aeea8c6`.
- Feature extractor: `deconv_features.py`, SHA-256 `b2cc85e294393c2d08382cf3a75d7ad9c9240a238dcd82f04f4e55c8f6bb35cb`.
- Synthetic tests: 5 passed; include scan-time unit regression.
- Frozen candidate and label files: hashes in `PLAN_V13b.code.sha256`.
- Raw features: `deconv_features_raw01_v13b.tsv`; evaluation: `DECONV13_EVALUATION_raw01_v13b.txt` under benchmark `astral-voxel-pilot-20261110/v2/v4/`.

## Next step
Close this deconvolution variant as a diagnostic failure under its predeclared threshold. Keep D11FWT unchanged. Do not promote to C++ or train on V13b. Resume the drawing-board analysis for a more entrapment-symmetric inference signal.
