# Raw 02 shape diagnostic: the raw 01 ID gain does NOT replicate, and calibration breaks further

Binary `f98c0dca2c`. Feature columns 1–16 are byte-identical to the raw 02 byion run (`CMP.txt`). Wall time 16:47, peak RSS 25.0 GB. Same evaluation and LDA replica as raw 01 (`shape_ablation.py`):

| arm (base = byion) | N at FDP_c 1% | plain q≤0.01 (FDP_c) | bulk-null q≤0.01 / 0.05 (FDP_c) |
|---|---|---|---|
| base | **1296** | 1236 (0.007) | 1115 (0.002) / 1663 (0.062) |
| +xcorr | 1172 | 1111 (0.007) | 1083 (0.007) / 1559 (0.061) |
| +library | 1282 | 1466 (**0.021**) | 1431 (**0.020**) / 2010 (**0.091**) |
| +signal | 1075 | 1300 (**0.019**) | 1272 (0.017) / 1723 (0.067) |
| +all (tool) | 1222 | 1406 (**0.020**) | 1353 (**0.021**) / 1877 (**0.073**) |

Reading:
- **No shape subgroup passes on both runs.** On raw 02 every subgroup loses entrapment-matched IDs at FDP_c 1% against the byion base.
- **The library-agreement subgroup inflates decoy-estimated IDs while the true error doubles.** On raw 02, q ≤ 0.01 has FDP_c 0.020–0.021 and bulk-null does not repair it. This replicates the CALIB section 8 finding: entrapment false targets are more library-like than decoys, so library-agreement scores mostly separate false targets from decoys and not true from false.
- **The raw 01 +22% was mostly this artefact plus run-specific variance.** The shape group stays implemented, opt-in and off.
- **Native design stays at search+boundary+byion.**
  - Raw 01: 1010–1029 against OSW 1151.
  - Raw 02: 1296–1298 against OSW 1287.
  - bulk-null q ≤ 0.01 FDP_c 0.006 / 0.002, and q ≤ 0.05 FDP_c 0.048–0.050 / 0.060–0.063.
