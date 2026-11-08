# Where the gap to OSW+Percolator comes from (raw 01, design diagnostic)

This compares OSW+Percolator and OpenDIA on the same library and the same raw. OSW gives 1151 at an entrapment FDP_c of 1%. We give 757 at matched FDP_c (`../matched_fdp.txt`).

The OSW-extracted XICs (`*_extracted.chrom.xic`, parquet, numpress) were decoded for all 588,890 precursors that map to our PQP. That gives 6 fragments, MS1 isotopes i0–i2, about 120 points per trace, and ±97 s around the library RT. They are the common ground for a model-free comparison.

DIA-NN's 1% IDs are used only as a reference: to define "confirmed right peak" and to bound what is learnable. They are never used as a pipeline signal.

## 1. Peak widths (`widths.txt`)

These are 1303 precursors confirmed by both DIA-NN and OSW, each with our candidate at the DIA-NN apex.

| | median | IQR |
|---|---|---|
| DIA-NN base width (RT.Start–RT.Stop) | **9.8 s** | 8.2–11.5 |
| DIA-NN FWHM | 3.6 s | 2.9–4.0 |
| OSW base width (left–right) | 13.1 s | 13.1–14.7 |
| our search group (start–end) | **45.8 s** | 19.7–76.9 |
| our boundary picker, base width | 15.1 s | 12.2–17.9 |
| our boundary picker, FWHM | 5.4 s | 5.1–6.0 |

- **Apexes agree.** Our apex minus DIA-NN's is 0.0 s at the median, and OSW's minus DIA-NN's is 0.4 s. The peaks are found. Only the group extent is wrong: just 20% of our group lies inside DIA-NN's span.
- **Why our groups are wide:** the streaming search (`OpenDIASearch.cpp::updateCandidate_`) closes a group only after a gap of more than 30 s between matched scans (`peak_gap_seconds = 30`), or when the valley split fires. Across all 2.8M groups the median is 67 s wide with 16 matched scans. That is 2.6 matched scans per 10 s, against 6.25 if every 1.6 s cycle matched. Groups are sparse streaks of single-fragment matches: 80% have a maximum transition support of 1 or 2. Their "width" is the span of scattered hits, not a chromatographic peak.
- **The boundary picker is already near DIA-NN/OSW** (15 s base, 5.4 s FWHM on 1.6 s cycles). Its widths are not the problem.

## 2. Our extraction matches OSW's (`extraction_cmp.txt`)

On the 512 anchors, our raw-scan traces against OSW's XIC interpolated onto our scans:

- Pearson correlation median 1.00 (q10 0.98);
- intensity ratio log2 −0.01;
- no scans where OSW is zero and we are not.

We are zero on about 12% of scans where OSW has small signal, because our calibration-search window is 15 ppm against OSW's 34 ppm. This catches the faint tail only.

**Limit of this check:** the anchors are bright, clean peaks, and the trace export uses the 15 ppm calibration window. The CNN and stage-2 traces use 0.02 Da (about 20–40 ppm), and the streaming search uses the calibrated 5 ppm tolerance. The weak, interfered precursors we miss were not compared. **On bright peaks, extraction is not the gap. For weak peaks this is untested.**

## 3. Which OSW IDs we miss (`ours_vs_xic.txt`)

These are 1462 OSW 1% IDs, 710 of them also accepted by us. The measures are taken on the OSW XIC inside OSW's boundaries.

| our recall of OSW IDs by quintile | Q1 | Q2 | Q3 | Q4 | Q5 |
|---|---|---|---|---|---|
| apex intensity | 0.04 | 0.28 | 0.51 | 0.72 | 0.88 |
| fragment correlation inside the OSW peak | 0.06 | 0.26 | 0.46 | 0.74 | 0.91 |
| share of trace signal inside the peak (low means background or interference) | 0.20 | 0.37 | 0.50 | 0.57 | 0.78 |

On the misses, the median fragment correlation is 0.64 against 0.93 for the ones we accept, and the in-peak signal share is 0.33 against 0.55. **OSW accepts weak, partly interfered peaks that we reject.**

## 4. Data treatment, using the same learner and the same candidates (`regimes2.txt`, `regimes3.txt`)

Candidates are the top-3 local maxima of the summed fragment trace within ±60 s of our calibrated library RT. Each regime changes one thing:

- **wide:** ±22 s boundaries;
- **picker:** valley or 5% boundaries;
- **baseline:** picker plus per-fragment rolling 20th-percentile subtraction.

The right peak is among the candidates for 2264 of 2731 confirmed precursors. The decoy-free check (AUC of confirmed right peaks against decoy groups) is high for co-elution (0.98) and shape (0.96–0.98) in every regime.

| supervised reference (bounded by 83% candidate coverage), q ≤ 0.01 / 0.05 (FDP_c at 1%) | wide 45 s | picker | baseline |
|---|---|---|---|
| MLP, 11 simple features | 961 / 1448 (0.025) | 1021 / 1432 (0.030) | 1013 / 1494 (0.020) |
| MLP, + 6 OSW-style xcorr / dot-product features | 806 / 1606 (0.038) | 432 / 564 | 844 / 1804 (0.024) |

- **Boundaries alone barely move the bound** (wide against picker: 961 against 1021 at 1%). With a good apex and a learner, a 45 s group still carries most of the signal. Narrow boundaries mainly improve width-sensitive features (width AUC 0.55 → 0.75, library correlation 0.80 → 0.86).
- **Baseline subtraction:** small and inconsistent at 1%, +4% at 5%. It is not the lever.
- **The 1% counts are unstable** across seeds and feature sets (for example 432 against 844), because the 1% cut sits in a sparse tail. The 5% column and the FDP_c are more informative.

## 5. Does a learned representation of the raw traces add information? (`cnn_oracle.txt`)

This uses the same candidates and folds. The input is 6 fragments + 3 MS1 isotopes × 33 scans (±26 s), max-normalised, plus log scale. The model is a supervised 1D-conv encoder plus an explicit 6×6 fragment correlation matrix, MLP head, 2 seeds.

| oracle | q ≤ 0.01 (confirmed right, FDP_c) | q ≤ 0.05 |
|---|---|---|
| (a) hand features, MLP | 808–1018 (FDP_c 0.014–0.017) | 1521–1600 |
| (b) trace CNN only | 464–658 (0.004–0.021) | 1201–1340 |
| (c) CNN + hand features | 838–1021 (0.012–0.018) | 1551–1694 |

With fully supervised labels, a small trace CNN does not beat 11 hand features, and adding it to them gains about 5% at 5%. The traces as extracted (6 fragments, ±26 s, max-normalised) do not hold much more separable information than simple summaries. That holds at least for this architecture and about 1900 positives per training fold.

## Reading: where the gap is

1. **Not peak finding or the apex, and extraction matches on bright peaks.** Our search has the right apex for 95% of confirmed precursors. The extraction on weak peaks is untested (section 2).
2. **Peak width is a symptom, not the main cause.** Our 45–67 s groups come from the 30 s gap rule over sparse single-fragment matches. Fixing them helps width-dependent features, but the oracle bound barely changes.
3. **Baseline subtraction is not a lever.** Interference matters a lot, though: the OSW IDs we miss are the low-correlation, low in-peak-share, low-intensity ones.
4. **The supervised reference on this harness reaches about 800–1000 at q ≤ 0.01, at FDP_c 1.4–3%.** At a true 1% that is below OSW's 1151, which OSW gets with decoy-only labels.

   This is **not an upper bound**. It is limited by the harness:
   - the candidates (top-3 summed-trace maxima, ±60 s) contain the right peak for only 83% of confirmed precursors, against 95% for our search;
   - the features are simple.

   The raw data are shared, so the gap must lie in what is computed from them: candidate generation, features, and score formation. Candidates that remain untested:
   - a per-candidate score from about 20 interference-aware subscores (S/N, library-weighted `var_xcorr_*`, isotope overlap, intensity share);
   - choosing the rank-1 peak group within a precursor by that score;
   - PyProphet- or Percolator-style iterative semi-supervised training.

   Which of these matters is tested next (`OSWLIKE.md`).
5. **Caveats:**
   - The supervised oracle uses DIA-NN-confirmed positives. Those are biased towards what DIA-NN can see and are a subset of true IDs.
   - OSW's 1151 includes precursors DIA-NN does not confirm, so the "confirmed-right" counts are not directly comparable with OSW's total.
   - The trace CNN was small and lightly trained (4 epochs). A larger model with more positives (cross-run) could do better. This test bounds the present setup, not representation learning in general.

## What this implies for the next step

The question "why do we lag OSW" is now narrower. It is about how many weak, partly interfered right peaks a score can lift above the decoy tail. The candidate levers, in order of evidence:

1. **Interference-aware evidence per fragment, not per group** (hypothesis). The misses have lower in-peak correlation and in-peak signal share. Their weakest fragment's share is 0.035 against 0.058. Whether one or two bad fragments explain this has not been shown.
2. **More labelled positives for the learned model.** In-run semi-supervised labels give about 700 positives; a cross-run model (raws 01 + 02) roughly doubles that.
3. **Fix the search grouping anyway** (cheap): close a group on a gap of 2–3 cycles rather than 30 s, so stage-2 features use real peaks. Expected gain is small.

The next measurable step is to compute OSW's own subscores for these exact candidates from the XICs and see which ones close the gap from about 1000 to 1151. That is an information question, not reuse of OSW's code.
