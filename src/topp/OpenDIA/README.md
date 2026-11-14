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

## Interference features (`OpenDIAInterference`)

Co-eluting precursors in one isolation window often share fragment m/z. A candidate that borrows a strong ion from its
neighbour then looks better than it is, and a real peptide whose fragments overlap a neighbour looks worse.

`OpenDIA::InterferenceResolver` links candidate peaks of different precursors that are in the same isolation window,
overlap in RT and share a fragment within the fitted tolerance. For every link, the claimant's abundance is estimated
from its own unshared fragments, and its elution profile gives how much of the shared signal it explains in the
candidate's window. That part is removed. Four features per candidate peak describe the result:

- `i_explained`: the explained share of the central signal;
- `i_libcos_resid`: the library cosine of the residual;
- `i_dlibcos`: the change in that cosine;
- `i_nclaim`: log(1 + number of claimants).

**Validation** (`gap-20261114`, PLAN_GAP arm I; Astral entrapment raws 01/02, 12-fragment DIALibGen library):
- Added to the stage-2 scorer, the Python prototype raised N at entrapment FDP 1 % from 1694 / 1628 to 1879 / 1921 (frozen seed).
- This C++ implementation gives +143..+185 (raw01) and +260..+367 (raw02) at 1 % over five training seeds, against +118..+185 and +199..+351 for the Python features.
- At 2 % the raw01 gain depends on the seed: +14..+131 for C++ and +39..+135 for Python. The raw02 gain is +205..+318.
- Entrapment FDP_c at q <= 0.01 stayed at or below 0.0071, and pooled decoy calibration (rho) was 0.97 [0.85, 1.11] at q <= 0.10.
- About a third of the raw02 gain comes from same-composition entrapment/parent pairs that only exist in the benchmark. A variant without those pairs still gains about +100 / +190.
- The C++ features equal the Python prototype given the same inputs (6 of 1.6 M windows differ by > 1e-5).
- The only intended input difference is the apex scan: the prototype's float binning puts it in either central bin, and here it always goes to the upper central bin.

**Cost:** 1:40 min and 3 GB per Astral raw; most of that is reading the mzML. The graph takes 9 s on 8 threads for 1.6 M candidate peaks and 54 M links, and 62 s on one thread. The output does not depend on the thread count.

## Tests

`ctest -R '^OpenDIA|OpenDIACalibrate'`. `OpenDIAInterference_test` checks four things:
- the analytic explained share;
- that windows without an RT overlap, a shared fragment or unshared claimant signal are not linked;
- equality with an all-pairs brute force on random data, and independence from the thread count;
- the extraction binning of the collector.
