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
- A 7.6K-parameter 1D segmentation model was trained on augmented anchors in about 1 min for 3 folds. On held-out anchors it beats the rule picker on the same binned windows, mainly at low SNR. IoU is 0.84 against 0.72 at high SNR and 0.66 against 0.43 at SNR < 1, with the rule run on scan bins only.
  - A first comparison ran the rule over empty bins and overstated the gap (0.84 against 0.52); it is withdrawn.
  - The windows are synthetic, so this is not evidence on real candidates.

**Not yet shown:** any gain in IDs.

**Stage-2 `boundary` group (raw 01, descriptive; plan written before the run; [`stage2-boundary/RESULT.md`](../gates/astral-raw01-peak-boundaries/stage2-boundary/RESULT.md); commit `5e23606c15`):**
- **What it is:** each window is re-picked on its scan bins, adding five features (log FWHM, log base width, co-elution and library correlation over the picked peak, apex agreement). No training is involved.
- **Counts:** `search+boundary` gives 360 / **746** / 1000 / 1441 at q ≤ 0.005 / 0.01 / 0.02 / 0.05, against 331 / 717 / 1015 / 1412 for the frozen `search`. That is +4.0% at 1%, below the pre-set +5%.
- **Checks:**
  - weights are stable in sign across folds;
  - null bands are symmetric;
  - the gain is on bright peaks, not weak ones;
  - CNN columns and HeLa defaults are byte-identical.
- **Decision:** no gate is proposed and the learned segmenter is not ported. The group stays opt-in. Better boundaries alone are not the main lever on this run.

**Next:** the boundary track is paused. The candidate levers from the miss analysis are the CNN's bias towards bright interference (the `-topaz_faithful` normalisation) and a non-linear stage-2 combiner. Each needs its own predeclared rule and an unused design run. A post-hoc raw 01 diagnostic ([LAG.md](../gates/astral-raw01-peak-boundaries/stage2-boundary/LAG.md)) ranks the combiner first. On the same search+boundary features, a 5-seed MLP ensemble reaches 935–1009 at q ≤ 0.01 against 745 for LDA, and 1534–1544 against 1440 at q ≤ 0.05, with R 0.98. Single-seed 1% counts are unstable (761–1131).
