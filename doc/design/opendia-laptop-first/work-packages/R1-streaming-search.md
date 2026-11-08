# R-WP1: Streaming sparse-grid candidate search

[Master plan](../README.md)

**Objective.** Generate high-recall candidate peak groups for the full library in one streaming pass over MS2 spectra, without materialising per-transition chromatograms.

**Delivered** (`feature/opendia-clean-rewrite`):
- `OpenDIA` reads the PQP library through the streaming importer, then packs spectra into m/z-sorted grid batches per isolation-window lane. Candidate blocks run in parallel (`8fdece8016`, `b47e3c8725`, `c63f33fe28`, `1ec6854adc`).
- Matching uses 1/K0 when present (`dcee451769`). `-im_array_as_vssc` covers unitless mzML IM arrays (`3944c4d552`).
- `-local_peak_split` splits peaks locally, and `-max_apexes_per_group` reports separated apexes (`cd5f50bec7`, `30bc59e4c1`). With n = 5, 13508 of the 19132 SVM IDs have a candidate apex inside the teacher bounds, against 6308 at n = 1.
- `-features_out` writes a feature table (`f4ecd4cb9e`).
- `-calibration` maps library coordinates to the run and extracts inside the calibrated RT, ppm and 1/K0 windows (`3671d01f94`). Candidates drop from 988,075 to 198,915 and extraction from 107 s to 40 s. Without the option the output is byte-identical.

**Open gates:**
1. Measure candidate recall at the calibrated windows: for teacher IDs with no in-window candidate, record why (RT, m/z, IM or support).
2. Recall must not regress on weak, missing-fragment and interfered precursors.
3. Native Thermo `.raw` through the same consumer.

## Peak-group boundaries, 8 November 2026 (raw 01 design data, diagnostic)

Evidence: [`gates/astral-raw01-peak-boundaries`](../gates/astral-raw01-peak-boundaries/RESULT.md), plus the miss analysis in [`MISSES.md`](../gates/astral-raw01-peak-boundaries/MISSES.md).

**Problem.**
- Against the 1462 same-library OpenSWATH+Percolator IDs on raw 01, the search finds a candidate at the right RT for 95% and picks the right one for 94%.
- 700 of these right peaks are not accepted at q ≤ 0.01.
- The search sets boundaries with a 30 s gap rule. Our median width is 3× OpenSWATH's on the same peaks (44 s against 13 s), which dilutes width, scans and every peak-integrated feature.

**Delivered (opt-in, defaults byte-identical, commit `e8cb5892b1`):**
- `OpenDIAPeakPicker`, a minimal consensus peak-group picker that takes about 240 ns per 6×32 window.
- `OpenDIACalibrate -out_anchor_peaks / -out_anchor_traces`, which picks each calibration anchor on its raw-scan traces to give teacher-free boundary labels.
- `augment_peak_labels.py`, which builds augmented training windows from those labels.

**Results:**
- All 512 raw 01 anchors give clean labels: FWHM 4.4 s, base width 14.7 s, right tail twice the left. These match OpenSWATH's widths, although OpenSWATH was not used to make them.
- A 7.6K-parameter 1D segmentation model was trained on augmented anchors in 59 s for 3 folds. On held-out anchors it beats the rule picker on the same binned windows at every SNR (IoU 0.84 against 0.52 at high SNR, 0.66 against 0.34 at SNR < 1).

**Not yet shown:** any gain in IDs.

**Next gate:** recompute peak features over the new [left, right], as one opt-in stage-2 group, with cross-fitted boundary models. Predeclare the rule and run it on an unused design run.
