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
