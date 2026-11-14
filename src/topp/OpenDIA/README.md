# OpenDIA

Identification-free DIA analysis on a predicted spectral library. Targets and decoys are treated identically; no labels
or external identifications are used.

## One command: `OpenDIAWorkflow`

```
OpenDIAWorkflow -in run1.raw run2.raw ... -lib library.pqp [-scoring_lib library12.pqp] -out_dir results -threads 8
```

Per run:
1. **Calibration**, identification-free: library RT to run RT map and tolerance, fragment m/z offset and fitted ppm
   tolerance, 1/K0 line, peak-boundary cap. A windowed second pass runs when the whole-run pass finds too few anchors.
2. **Candidate search** with the calibrated library: native-scan XIC peak picking, up to 5 peaks per precursor.
3. **Population**: every candidate of each target-decoy pair whose best candidate has at least 3 supporting fragments.
4. **One more pass over the spectra** collects the scoring evidence: 6-fragment traces and apex b/y ions, and the
   12-fragment grid of every window. A grid is held only while its window is open.
5. **24 features**:
   - search 4;
   - boundary 5;
   - b/y ions 2;
   - 12-fragment evidence 6;
   - fragment deficit 3;
   - shared-signal interference 4.
6. **Scoring**: `OpenDIA::CrossFitScorer`, three folds by peptide family. Seeds are the calibration anchors' windows
   and decoy windows. The scorer is a small network with a decoy-whitened readout and two self-training rounds. Then
   comes one window per precursor, decoy-tail calibrated scores, and pooled paired target-decoy q-values.

Across runs, a precursor's global q-value comes from its best run score, by paired target-decoy competition.

**Outputs** (in `-out_dir`):
- `<run>.precursors.tsv`: one row per precursor, with its window, score and run q-value;
- `<run>.calibration.json`;
- `precursors.tsv`: target precursors at run and global q ≤ `-q`, all runs;
- with `-features_out`, `<run>.features.tsv`: every scored window with its 24 features.

**Raw files** need a .NET 8 runtime (set `DOTNET_ROOT` if it is not in a standard location); mzML works without it.

**Cost** on one Astral 15 min raw file (8 threads): about 8-10 min and 9.6 GB peak.
- Search: 30-40 s.
- Evidence pass: 5-6 min. It holds about 90,000 open 12-fragment grids and the 6-fragment trace tensor of 1.6 M windows.
- Scorer: 1-2 min.
- Calibration adds 1-3 passes when it is not reused.

## Single steps

The step tools stay for inspection and experiments.

| Tool | Step |
|---|---|
| `OpenDIACalibrate` | Calibration only (`cal.json`, `-out_anchors`) |
| `OpenDIA -calibration cal.json -features_out` | Candidate search and its feature table |
| `OpenDIAInterference` | Interference features of given windows |
| `OpenDIACNN` (optional, `-DOPENDIA_WITH_TORCH=ON`) | The earlier CNN + LDA scorer. It needs LibTorch and is not part of `OpenDIAWorkflow`. |

## Interference features

Co-eluting precursors in one isolation window often share fragment m/z. A candidate that borrows a strong ion from its
neighbour then looks better than it is, and a real peptide whose fragments overlap a neighbour looks worse.

`OpenDIA::InterferenceResolver` links candidate peaks of different precursors that:
- are in the same isolation window;
- overlap in RT;
- share a fragment within the fitted tolerance.

For every link, the claimant's abundance is estimated from its own unshared fragments. Its elution profile gives how
much of the shared signal it explains in the candidate's window, and that part is removed. Four features describe the
result:
- `i_explained`: the explained share of the central signal;
- `i_libcos_resid`: the library cosine of the residual;
- `i_dlibcos`: the change in that cosine;
- `i_nclaim`: log(1 + claimants).

**Validation** (`gap-20261114`, PLAN_GAP arm I; Astral entrapment raws 01/02, 12-fragment DIALibGen library):
- **Python prototype:** added to the stage-2 scorer, it raised N at entrapment FDP 1 % from 1694 / 1628 to 1879 / 1921 (frozen seed).
- **C++ implementation:**
  - +143..+185 (raw01) and +260..+367 (raw02) at 1 % over five training seeds, against +118..+185 and +199..+351 for the Python features.
  - At 2 % the raw01 gain depends on the seed: +14..+131.
  - Entrapment FDP_c at q ≤ 0.01 stayed ≤ 0.0071.
- **Benchmark-specific share:** about a third of the raw02 gain comes from same-composition entrapment/parent pairs that only exist in the benchmark.

## Tests

Run `ctest -R '^OpenDIA|OpenDIACalibrate'`.
- **`OpenDIAScorer_test`:**
  - gradients against finite differences;
  - the readout in closed form;
  - paired competition on a hand example;
  - tail calibration monotonicity;
  - thread independence and recovery of a synthetic population.
- **`OpenDIAFragmentGrid_test`:** the 12-fragment and deficit formulas against numpy reference values, and the streaming grid collector against a brute-force fill.
- **`OpenDIAInterference_test`:** the explained share against an all-pairs brute force.
