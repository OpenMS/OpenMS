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

## 6. Fix 2 applied: RT window at the gradient ends (OpenMS commit 8d2771bcec)
- **Cause.** The anchors come from the windowed second pass. Its windows are centred on the stage-1 map, which at library RT < 0.05 points to 450 to 490 s. So no panel precursor eluting at 190 to 340 s was collected, and the final map extrapolated there (CAL_EARLY_raw01.txt, STAGE1_EARLY_raw01.txt).
- **Part A, collection.** Panel members whose library RT is outside the stage-1 anchor support are collected and searched to the run edge (89336 targets). On its own (rtedge/) this found 5 correct early anchors, but the 20-bin map rejected them as outliers.
- **Part B, map.** Stage 4 now fits an L1 isotonic map (one knot per PAVA block), extrapolated with Theil-Sen end slopes over sqrt(n) end anchors. Sigma is the 10-fold cross-validated residual. Outliers are rejected on that held-out residual. The cal.json rt block gains map, low_slope_s, high_slope_s; old files still read.
  - Rejecting on the in-sample residual instead was tried first (rtmono_insample/). It kept a lone late anchor (library RT 0.869, 1321 s, not in DIA-NN's report) that bent the late end. Held-out rejection drops it (rtmono_reject.py, REJECT_VARIANTS_raw01.txt).
- **raw01 result** (rtmono/RTEDGE_EVAL_raw01.txt):

  | | Before (fitppm) | After (rtmono) |
  |---|---|---|
  | Used anchors | 527 | 529 |
  | Anchor run RT range | 409 to 1182 s | 188 to 1182 s |
  | RT sigma / window half width | 21.0 / 84.1 s | 21.1 / 84.4 s |
  | DIA-NN q <= 1% apexes outside our window | 62 (3 entrapments) | 30 (2 entrapments) |
  | Of those, library RT < 0.06 | 48 | 14 |
  | Median abs(DIA-NN apex - map), library RT 0.02 to 0.04 | 48.9 s | 7.3 s |
  | Apex +-10 s inside window, q <= 1% | 96.3% | 97.9% |
  | Missed DIA-NN points / extra points, q <= 1% | 3.9% / 0 | 4.0% / 0 (over 1.7% more compared points) |
  | Calibration / search | 3:21, 16.8 GB / 2:42, 4.4 GB | 3:51, 20.3 GB / 3:10, 4.4 GB |

- **What remains.** 30 apexes outside. Library RT < 0.02 still has a median map error of 55 s (3 anchors there). The late end (library RT > 0.8, 8 misses) has 12 anchors. Calibration memory went up 3.5 GB from the run-edge collection.

## 7. Peak parity: apex, boundaries and peak choice (raw01, rtmono run)
Scripts: peak_parity.py, boundary_rule.py, diann_peak_rule.py, our_boundary_rules.py, peak_choice.py, picker_port.py. Outputs in rtmono/ (PEAK_PARITY, BOUNDARY_RULE, DIANN_PEAK_RULE, OUR_BOUNDARY_RULES, PEAK_CHOICE, PICKER_VARIANTS _raw01.txt). Scan cycle 1.646 s. DIA-NN's RT, RT.Start and RT.Stop sit on our scan RTs (|d| <= 0.0001 s), so all comparisons are in whole scans.

**Apex and peak choice, DIA-NN q <= 1% (2357):**

| | Count |
|---|---|
| DIA-NN apex outside our RT window | 29 (2 entrapments) |
| Inside the window, but in none of our 5 reported peaks | 30 (2 entrapments) |
| Contained in one of our 5 peaks | 2298 (97.5%) |
| ... and that peak is our rank-0 | 2133 (92.8% of contained) |
| Apex of the matched peak: same scan / 1 scan / 2 / 3-6 / >6 | 1655 / 542 / 55 / 49 / 56 |

- The apex agrees: median offset 0, 10-90% range -1 to +1 scan. DIA-NN's own apex is the raw max of its fragment sum in 72% of cases and the [1 2 1]-smoothed max in 76% (diann_peak_rule.py), so a 1-scan spread is inherent to its rule, not to our data.
- The FWHM agrees: ours / DIA-NN median 0.88.

**Boundaries (peaks containing DIA-NN's apex, n = 2298):**
- Start: exact 52%, median offset 0.
- Stop: exact 13%, median +3 scans, p90 +7. Our width median is 9 scans, DIA-NN's 6 (max 12; it never goes beyond 6 scans per side, matching "Scan window radius set to 6" in its log). Our peak contains DIA-NN's [start, stop] in 63% of cases.
- Cause: our right boundary walks down the elution tail to 5% of the apex. DIA-NN's stop sits on average at 6-8% of the apex, and the trace keeps falling after it (median next scan / stop scan 0.72), so DIA-NN cuts the tail early. 19% of our long ends also run through a merged shoulder (valley above 0.5 of the lower maximum).
- Rules on our traces (our_boundary_rules.py and a FWHM cap check): capping each side at 1.25x our own FWHM gives stop exact 38%, within 1 scan 78%, median 0, IoU 0.83 (now 0.67). Stopping at the first smoothed minimum gives stop exact 21%. A 10% floor plus a 6-scan cap gives stop exact 32%, IoU 0.83.

**Peak choice: a degenerate co-elution score.** 20% of rank-0 peaks (targets and decoys alike) have co-elution exactly 2.000, all with transition_support 2 and median matched_scans 1. Two fragments that are nonzero on the same single scan correlate with r = 1 over any window, so a one-scan coincidence ties a 2-fragment perfect elution. Of the 165 q <= 1% precursors where DIA-NN's peak is one of our 5 but not rank-0, 51 have an integer rank-0 co-elution. Only 1 of the 2135 correct rank-0 picks does.
- The Python port of pickXicPeaks reproduces our reported peaks for 2314 of 2357 precursors. 43 of the 45 differences are ties at 2.000 ordered differently, so the port is faithful.
- Requiring a fragment to have signal on >= 2 scans of the correlation window before it enters co-elution (profile2): rank-0 = DIA-NN's peak goes 2133 -> 2189 (+56), contained 2298 -> 2304, entrapments among the rank-0 matches 15 -> 15. Boundaries unchanged.
- Cutting boundaries at the first smoothed minimum loses 37 contained peaks (the cut peak no longer contains DIA-NN's apex). It is worse than a FWHM-scaled cap.

**Not changed in code yet.** Both candidate fixes (profile2 co-elution; a FWHM-scaled boundary cap) change the candidate population that the scorer sees, so they need a decision and a frozen plan before a run.

### 7a. Does a width-scaled boundary cap generalise? (cap_generalisation.py, rtmono/CAP_GENERALISATION_raw01.txt)
- raw01 covers only a narrow range of widths: our FWHM tertile medians are 1.4 / 1.9 / 2.6 scans. Wider peaks (longer gradients, slower cycles) are tested only by simulation.
- **Matching DIA-NN, by width tertile:** the k that matches best changes with width (sym: 1.5 / 1.25 / 1.0-1.25 for narrow / mid / wide). DIA-NN's width / our FWHM is 3.7 / 2.9 / 2.5. So DIA-NN is not scale-free itself: it has a minimum width at the sampling limit plus a fixed 6-scan radius set from its run-level peak width. Any single k matches DIA-NN less well at the extremes of this run.
- **By tailing tertile (hr/hl 0.56 / 1.15 / 2.4):** the symmetric cap fits tailing peaks at k 1.5, but it cuts the tail through the left width. The per-side cap (left k x hl, right k x hr) is the shape-aware form. Its best k is 2.0-3.0 across all strata.
- **Simulation** (Gaussian and EMG, FWHM 1.5 to 20 scans, counting noise, cap from the noisy FWHM estimate). Share of the true area inside the cap:
  - Gaussian: 99.6-100% at every width.
  - Tailing tau/sigma 1: 97.7-99.3% (sym1.25), 98.0-99.3% (side2.5).
  - Tailing tau/sigma 2: 92.5-96.2% (sym1.25), 94.6-97.2% (side2.5).
  - The width dependence is small (a slight loss at wide tailing peaks). The tailing dependence dominates.
- **Conclusion:** a cap in units of the peak's own half widths is scale-free in width. It loses only the far tail of strongly tailing peaks, and less so per side. The k value itself (1.25 FWHM, 2.5 half widths) was chosen on raw01 against DIA-NN. It is a dimensionless shape constant like the existing 5% floor, but it is not derived from the data.
