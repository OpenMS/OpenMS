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
