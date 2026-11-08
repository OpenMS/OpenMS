# Raw 01 miss analysis (descriptive, diagnostic only)

Scripts: `miss01.py`, `band01.py`. Outputs: `miss01.txt`, `band01.txt`.
- Nothing was fitted or selected, and entrapment labels were not joined.
- The references are used only to locate misses:
  - DIA-NN 1% exact library matches, 2754 precursors;
  - same-library OpenSWATH+Percolator 1%, 1462 precursors.
- Our result is the frozen grid2d CNN + stage-2 search set: 717 targets at q ≤ 0.01.

## Funnel

| step | vs OSW 1462 | vs DIA-NN 2754 |
|---|---|---|
| a candidate within 10 s of reference RT | 1394 (95%) | 2449 (89%) |
| our winner at reference RT | 1377 (94%) | 2292 (83%) |
| accepted at q ≤ 0.01 | 680 (47%) | 702 (25%) |
| accepted at q ≤ 0.05 | 1106 (76%) | 1285 (47%) |

Losses against OSW:
- 68 have no candidate at the reference RT;
- 17 have a right candidate that lost the winner slot;
- **700 have the right winner but are not accepted.**

## Who is missed
- **Near misses:** the median q of the 700 is 0.036, and 25% sit at q ≤ 0.018.
- **Weak signals:** recall of OSW IDs by OSW intensity tertile is 13% / 49% / 77%. By OSW d_score quartile it is 12% / 30% / 60% / 84%.
- **Charge is not the cause:** charge mix matches OSW.

## Why (medians; AUC = missed right-peak vs top 1% decoy winners)

| feature | ours accepted | OSW right peak, missed | top 1% decoys | AUC |
|---|---|---|---|---|
| cnn_z | 10.5 | 6.2 | 5.0 | 0.64 |
| transition support | 6 | 5 | 4 | 0.84 |
| library correlation | 0.93 | 0.85 | 0.17 | **0.86** |
| fragment co-elution | 0.88 | 0.54 | 0.14 | 0.82 |
| MS1 ppm | 0.6 | 1.6 | 7.2 | 0.25 |

- **The CNN barely separates these** (AUC 0.64).
- **Library correlation and co-elution separate well** but are not in the frozen search set.
  - In the q 0.01 to 0.05 band (698 targets, 66 decoys), library correlation has within-band AUC 0.77.
  - 237 of 698 band targets and 3 of 66 band decoys have co-elution > 0.6 and library correlation > 0.6.
  - The earlier linear single-feature add (ablation) did not convert this into IDs: library correlation gave 705 against 717.
- **Peak boundaries are too wide:**
  - our median width on OSW-matched right peaks is 44 s, against 13 s for OSW (ratio 3.1);
  - stage-1 FWHM is about 3.3 s;
  - matched scans cover about 0.39 per second of width.
  - Width, scans and window-integrated features are therefore diluted.
- **What band decoys look like:** brighter (log apex 4.9 vs 4.4), narrower (20 s vs 36 s), further from library RT (32 s vs 16 s), and at lower m/z (432 vs 546).
  - This is bright co-located interference.
  - The m/z asymmetry is a warning: m/z must not become a feature.

## What is missing against OpenSWATH
- **Peak-group boundaries** from a real picker.
- **About 40 subscores over those boundaries:**
  - library dot product and spectral angle;
  - cross-correlation co-elution and shape, weighted by library intensity;
  - per-transition S/N;
  - MS1 isotope overlap;
  - normalised RT.
- **A second stage that combines them.**
- Our 598-seed CNN plus 5 search features is not enough. Six transitions are not the limit, since OSW uses the same 6.
