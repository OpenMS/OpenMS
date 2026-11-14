# OpenDIA

Identification-free DIA search on a predicted spectral library. The tools run in this order. Every tool streams the
raw file once and treats targets and decoys identically. No labels or external identifications are used.

| Step | Tool | Output |
|---|---|---|
| 1 | `OpenDIACalibrate` | `cal.json`: monotone library-RT to run-RT map and RT tolerance, fragment m/z offset and fitted tolerance, 1/K0 line, peak-boundary cap; `-out_anchors` |
| 2 | `OpenDIA -calibration cal.json -features_out` | Candidate peaks per precursor (native-scan XIC peak picking, up to `-max_peaks`) |
| 3 | `prepare_trace_windows.py --anchors` | Scoring windows and seed windows from the calibration anchors |
| 4 | `OpenDIACNN -mode crossfit -stage2 -stage2_feature_set search+boundary+byion -paired_tdc -pooled_tdc -rt_calibration cal.json -apply_calibration_coordinates` | Cross-fitted CNN logit plus a semi-supervised LDA over the CNN logit and peak-group features; pooled paired target-decoy q-values |

`OpenDIACNN` needs a C++11-ABI CPU LibTorch: configure with `-DOPENDIA_WITH_TORCH=ON -DTorch_DIR=...`. The other tools do not link LibTorch.

## Tests

`ctest -R '^OpenDIA|OpenDIACalibrate'`.
