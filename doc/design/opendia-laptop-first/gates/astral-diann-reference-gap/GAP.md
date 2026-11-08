# OpenDIA vs ProteoBench DIA-NN 2.3.1 on the Astral entrapment runs (descriptive)

DIA-NN is used only as a reference here. Nothing in the pipeline was changed or fitted. Raw 02 is spent, and raw 03 was only read through DIA-NN's own output.

## Libraries differ about 30-fold
- **DIA-NN:** a full predicted library of 8.65M precursors (2.69M peptides, about half entrapment). It finds about 82K precursors per run at 1% (82360 / 82498 / 81685, from `param_0..txt`).
  - ProteoBench's `result_performance.csv` keeps 18045 / 53239 / 22377 per run after its own filtering.
- **Ours:** 100,000 peptides randomly sampled from the same entrapment FASTA (49,683 entrapment), giving 300,223 target precursors.
  - Charge 1–5 with 6 transitions each; 124K of the precursors are charge 4–5.
- **Overlap:** DIA-NN 1% precursors that exactly match one of our precursors (sequence, charge and m/z including modifications) number **2754 / 2745 / 2757** per run (11 / 14 / 8 entrapment). Of these, 74% are charge 2 and only about 50 are charge 4.
- **Ceiling:** about 2750 per run is the realistic one for this library. The 90K figure is not reachable with this library.

## Funnel on raw 02 (2745 DIA-NN-confirmed precursors in our library)

| step | precursors |
|---|---|
| any search candidate | 2734 |
| a candidate within 10 s of DIA-NN RT | 2449 |
| stage-2 winner at that RT | 2215 |
| accepted at q ≤ 0.01 | **500** (518 accepted in total, 506 DIA-NN-confirmed) |
| accepted at q ≤ 0.05 | 876 |
| accepted at q ≤ 0.1 | 1014 |

- **Recovered:** 18% of the ceiling. On HeLa we reached about 70% of OpenSWATH+PyProphet.
- **Main loss is scoring:** 1715 precursors have the right peak chosen but are not accepted. Their medians are close to the top 1% of decoy winners on most features.

  | | cnn_z | co-elution | library correlation | MS1 ppm |
  |---|---|---|---|---|
  | right peak, not accepted (1715) | 6.3 | 0.46 | 0.73 | 2.1 |
  | top 1% of decoy winners | 6.0 | 0.23 | 0.37 | 7.2 |
  | accepted | 10.9 | 0.89 | 0.92 | 0.5 |

  - Library correlation and MS1 ppm still separate, but co-elution and cnn_z largely overlap.
- **Smaller losses:**
  - search: 285 precursors with no candidate within 10 s;
  - winner choice: 234 precursors where a wrong peak won.
- **RT:** DIA-NN RT against library RT has a robust residual of 28.6 s (q90 61.7 s), so the library RT itself is coarse. Our calibration σ of 20.7 s is consistent with that.

## Same-library baseline: OpenSWATH + Percolator (`100K_sampled_fasta/opendia_results`, older OpenDIA workflow, all 3 runs)

| run | accepted at 1% | entrapment | FDP_c | DIA-NN-confirmed |
|---|---|---|---|---|
| 01 | 1441 | 23 | 0.032 | 1361 |
| 02 | 1470 | 25 | 0.034 | 1395 |
| 03 | 1396 | 22 | 0.032 | 1317 |

- **Calibration:** the baseline is anti-conservative too, at about 3.3% entrapment FDP for a nominal 1%.
- **Matched entrapment FDP (about 0.03, raw 02):** stage 2 accepts about 776 against the baseline's 1470. Our pipeline recovers about half of it.

## Reading
- **The ~90K figure does not apply to this library.** About 2750 precursors per run is the DIA-NN-confirmed ceiling. The reference workflow on the same library reaches about 1400 confirmed.
- **The gap to close:** 500 confirmed against about 1400 for OpenSWATH+Percolator, out of a ceiling of 2750.
- **Where the gap is:** scoring power on peaks we already find (1715 right-peak precursors not accepted). Search and winner choice matter less.
- **Plausible causes, not tested:**
  - the CNN is trained on only 598 positive seeds;
  - stage 2 has 8 features, against the dozens of OpenSWATH subscores;
  - there are only 6 transitions per precursor.
