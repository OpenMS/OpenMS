# Scoring ablation on raw 01 (design run, target/decoy labels only)

All numbers here are pooled paired q from an offline replica of `OpenDIACNN -stage2` (`lda.py`). The replica reproduces the tool exactly: 333 / 525 / 770 at q ≤ 0.01 / 0.02 / 0.05.
- Same cross-fitted semi-supervised ridge LDA, with per-fold decoy standardisation.
- No settings were tuned, and entrapment labels were not joined.
- DIA-NN agreement (full-library ProteoBench run) appears only as a diagnostic column. It is never used to fit or select.

## What the scorer sees today
- **CNN input:** 6 MS2 fragment channels × 32 RT bins (±15 s) × 8 IM bins, plus library intensities and an IM flag.
  - No MS1 trace and no isotopes.
  - On Astral (no IM) every point lands in the middle IM bin, so 7/8 of the tensor is zero.
  - Model: 2 conv layers (8 and 12 channels) and a 24-unit head.
  - Trained on 598 positive seeds per run.
- **Stage 2:** MS1 monoisotope trace correlation, MS1 ppm, and presence of M+1 in a ±4 s core window. Nothing beyond M+1.
- **For comparison, OpenSWATH:**
  - extracts MS1 M+0..M+3 (`OpenSwathWorkflow -ms1_isotopes 3`);
  - scores isotope correlation and overlap with `dia_nr_isotopes 4`, at MS1 and MS2;
  - has about 40 subscores in total, including `var_norm_rt_score` (library ΔRT).
- **Not used:** the search's own peak features (transition support, matched scans, peak width, library ΔRT) were in `search.features.tsv` but were never given to the scorer.

## Feature groups added to the frozen stage-2 LDA (counts at q ≤ 0.01 / 0.02 / 0.05)

| features | counts | share confirmed by DIA-NN (same RT ±10 s) at 0.01 / 0.05 |
|---|---|---|
| cnn_z only (frozen CNN-only design) | 321 / 462 / 634 | n/a |
| stage 2, frozen (8 features) | 333 / 525 / 770 | 0.976 / 0.906 |
| + library ΔRT (abs, squared/σ) | 363 / 537 / 784 | 0.970 / 0.913 |
| + intensity (log apex, log area) | 346 / 553 / 759 | n/a |
| + charge one-hot | 327 / 526 / 771 | n/a |
| + search shape (support, log scans, log width) | 623 / 913 / 1316 | 0.979 / 0.904 |
| + search shape + ΔRT | 646 / 1024 / 1343 | 0.977 / 0.901 |
| **cnn_z + search shape + ΔRT (no extraction features)** | **718 / 1014 / 1411** | n/a |

## Single extraction features added to cnn_z + search shape

Baseline (cnn_z + search shape) is 655 / 875 / 1331.

| added feature | q ≤ 0.01 / 0.02 / 0.05 |
|---|---|
| abs ΔRT | 717 / 1014 / 1412 |
| library correlation | 705 / 918 / 1366 |
| ms2 ppm | 668 / 828 / 1292 |
| ms1 ppm | 656 / 901 / 1323 |
| ms1 isotope | 655 / 835 / 1331 |
| fragment co-elution | 626 / 918 / 1383 |
| ms1 co-elution | 626 / 990 / 1341 |
| fragment fraction | 622 / 803 / 1211 |

## Symmetry check
- **Question:** a feature could exploit target/decoy asymmetry rather than real signal. Among CNN-z bands that are mostly false (z ≤ 6), do the search features differ between targets and decoys?
- **Result:** they barely differ. Mean transition support is 1.52 vs 1.48 at z ≤ 0 and 2.91 vs 2.78 at z 4 to 6, and widths and scans are similar.
- **Diagnostic agreement:** DIA-NN agreement of accepted IDs stays at about 0.977 at q ≤ 0.01 and about 0.90 at q ≤ 0.05, the same as the frozen design.
- **Limit:** this is not entrapment evidence.

## Reading
1. **Largest lever:** the search's own peak-group evidence, chiefly transition support. With ΔRT it doubles the q ≤ 0.01 count, from 333 to 718. That is still about half of the same-library OpenSWATH+Percolator baseline (1441 on raw 01).
2. **Library ΔRT:** a small, consistent gain (+30 to +60).
3. **Current hand extraction features:** they add little once search features are in, and several reduce the count. Their information overlaps with transition support. They are crude: one MS1 isotope and Pearson correlations on 1 s bins.
4. **Not yet tested, because new extraction is needed:**
   - MS1 M+0..M+3 traces;
   - isotope correlation and overlap;
   - an MS1 channel inside the CNN;
   - a CNN without the empty IM axis;
   - more CNN positives through iteration on stage-2 labels.
