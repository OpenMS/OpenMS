# Extraction parity, OpenDIA xic vs DIA-NN 2.2.0, raw01 (2026-11-12)

DIA-NN 2.2.0 rerun on our library with `--xic 60 --export-quant`. Otherwise it is the same command as diann-ourlib-20261112.
- It reproduces 2357 IDs at q <= 1%. Wall time 11:03.
- It exports per-fragment and MS1 traces for all 7553 reported precursors (q <= 0.5), +-60 s around each apex.

OpenDIA ran with `-peak_picking xic` plus the new `-xic_traces_out` dump (commit f8f4742f40) for the same 7553 precursors.
- Wall time 3:17, max RSS 4.44 GB.
- Feature output is byte-identical to xic-peaks-20261112.

Scripts: `parity.py`, `parity2.py`, `miss_ppm.py`, `signal_ppm.py`, `ppm_vs_intensity.py`. Read-outs are in the `*_raw01.txt` files.

## 1. What already matches
- **Scans.** Every DIA-NN point falls on one of our native scans at the same RT (|dRT| <= 0.0001 s). DIA-NN exports one point per cycle (1.64 s), like us.
- **Intensities.** Where both tools have a value, it is the same raw peak: ratio 1.0000 for 100% of points.
- **Fragments.** Every exported DIA-NN fragment is one of our library fragments (34979 of 34979 paired).
  - DIA-NN exports only 4.6 of the 6 library fragments.
  - It drops every b1/b2/y1/y2 ion. Precursors with 0, 1 or 2 such ions export exactly 6, 4 or 3 fragments (FRAGSET_raw01.txt), so the rule is deterministic.
  - We extract all 6, which is a superset.

## 2. Difference 1: fragment m/z tolerance (the main cause)
- Ours: calibrated +-5.017 ppm (5 x 1.00 ppm sigma). DIA-NN: "Optimised mass accuracy: 8 ppm".
- At the original 5.017 ppm, on paired traces and scans:
  - 16% of DIA-NN's nonzero points for q <= 1% IDs are zero in our traces, and 13% within 5 s of the DIA-NN apex.
  - Per-trace Pearson r against DIA-NN: median 0.987, but 27% of traces are below 0.9.
- Every missed point is a raw peak at 5.02 to 8.1 ppm from our calibrated query (3000 of 3000 sampled, MISS_PPM_raw01.txt).
- **Probe:** the same calibration with only the half width set to 8.0 ppm (probe8/). Not a fix, a hypothesis test.
  - Missed DIA-NN points drop to 0.13% for q <= 1% IDs (0.10% within 5 s of the apex).
  - Our extra points are 0.44%. Both residuals sit at |ppm| 7.84 to 7.98, which is rounding at the edge of the window.
  - Per-trace r is 1.0000 at the median and above 0.99 for 91.5% of traces.
  - Cost: 3:34 wall, 4.43 GB.
- **Conclusion:** with the same tolerance, the grid extractor returns DIA-NN's fragment data point for point.
- **Why our sigma is too small** (SIGMA_ESTIMATES_raw01.txt, PPM_VS_INTENSITY_raw01.txt):
  - Calibration estimates sigma from an intensity-weighted mean m/z over +-15 ppm at the apex scan of each anchor only.
  - Single centroids scatter more, and the scatter depends on intensity. Robust sigma is 2.9 ppm in the weakest quintile and 0.7 ppm in the strongest.
  - Near the apex of confident IDs, 6.2% of real signal points lie beyond 5.017 ppm: 15% in the weakest decile, 0.2% in the strongest.

## 3. Difference 2: RT window outside the calibration anchors
- Our window contains the DIA-NN apex +-10 s for 96.3% of q <= 1% IDs (89% overall). 63 of the 2357 IDs have their apex outside our window.
- Cause: the calibration has no anchor eluting before 409 s.
  - Anchor library RT spans 0.026 to 0.864. 5.4% of library targets lie below that range and 16.7% above it.
  - Below the first knot the map extrapolates the end slope, putting library RT 0.01 at about 424 s. DIA-NN's IDs there elute at 190 to 360 s (RT_EARLY_raw01.txt).
  - 46 DIA-NN q <= 1% IDs lie below the first anchor (apex 191 to 418 s) and 6 lie above the last one.
- Our anchors themselves agree with DIA-NN: median |dRT| 0.0 s, max 8.2 s (ANCHORS_VS_DIANN_raw01.txt). The map is fine inside the anchor range. It is unsupported outside it.

## 4. Not compared yet
- MS1 traces (DIA-NN exports MS1; we do not extract MS1 in the search).
- Peak boundaries after the fix.

## 5. Fix 1 applied: fitted fragment tolerance (OpenMS commit 92cc9a2469)
- **Calibration change.** It now models the error of the strongest centroid at the apex spectrum and its two neighbours (8937 errors from 527 anchors).
  - The model is a two-Gaussian scale mixture plus uniform background over the +-15 ppm search window.
  - The half width is where signal density falls to background density. There is no sigma multiple.
- **raw01 result.** Offset -3.27 ppm. Half width 7.13 ppm, previously 5.02. Anchors and RT map are byte-identical. Calibration 3:21 / 16.8 GB, search 2:42 / 4.42 GB (fitppm/).
- **Parity against DIA-NN, q <= 1% IDs** (fitppm/TRACE_PARITY_raw01.txt):
  - Missed DIA-NN points: 16.4% before, 3.9% now. Within 5 s of the apex: 12.9% before, 2.7% now.
  - Extra points: 0.
  - Per-trace r > 0.99: 44% of traces before, 76% now.
- **What remains is the 7.13 to 8 ppm band.** DIA-NN uses 8 ppm, and the crossing depends on the model (TAIL_MODELS_raw01.txt, GRENANDER_raw01.txt):

  | Model | Crossing (ppm) | Note |
  |---|---|---|
  | 1 Gaussian | 4.1 | |
  | 2 Gaussians | 7.13 | lowest BIC |
  | 3 Gaussians | 8.12 | BIC +2, a tie |
  | Student-t | 9.05 | |
  | shape-free monotone density | 7.48 | bootstrap 7.2 to 7.5 |

  The 2-Gaussian model underfits the observed tail at 6 to 8 ppm (FIT_TAIL_CHECK_raw01.txt).

## 6. Not fixed yet: RT window at the gradient ends
- The anchors come from the windowed second pass. Its windows are centred on the stage-1 map, which at library RT < 0.05 points to 450 to 490 s.
- So no panel precursor eluting at 190 to 340 s is ever collected. The final map then extrapolates there (CAL_EARLY_raw01.txt, STAGE1_EARLY_raw01.txt).
- **Search-side option, simulated on DIA-NN q <= 1% apexes** (RT_RULES_raw01.txt): outside the anchor library-RT range, extend the window to the run edge.

  | Rule | DIA-NN apexes outside our window | Entrapments among them | Mean window width |
  |---|---|---|---|
  | Current | 61 | 3 | 1.0x |
  | Extend outside anchor range | 26 | 1 | 1.49x (22% of precursors widened) |
  | Extend outside knot range | 8 | 1 | 1.96x |

- **Calibration-side option:** let the windowed pass also collect panel members beyond the stage-1 map support up to the run edge, so anchors exist there.
