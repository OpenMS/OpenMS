# PLAN_CALREFINE: one stage-3 pass around the final map

Written 15 November 2026, before implementation. It follows PLAN_CALSPREAD, which failed offline.

## Problem
The stage-3 window is centred on an earlier map: the stage-2 EM map, or the stage-1 map in the windowed pass. Stage 4
fits the RT spread with:
- target residuals to the final map;
- a background taken in window coordinates.

When the final map differs from the window centre (median 11.4 s on arm E raw01), the two frames disagree and the
spread is inflated. Measured against DIA-NN 1 % IDs:
- arm E raw01: 37.7 s against 27.2 s;
- arm E raw02: 30.9 s against 25.4 s.

## Design (one variant, no tuning, no new parameter)
After stages 3-4 succeed, run stages 3-4 once more on the traces already collected, with the window centred on the
final map of the first pass. Everything else stays the same:
- the half window;
- the support bounds;
- the anchor rule;
- the spread fit.

The window centre then equals the map the residuals are taken to, so the existing fit's frames match. If the second pass
throws (too few anchors), the first result is kept. No extra raw pass is needed. In the windowed pass the collected
traces extend 30 s + 2 FWHM beyond the first window, which covers the shift on these runs; members whose recentred window
leaves the collected range are searched over what was collected.

The JSON records `stage3.refined_on_final_map`.

## Evaluation (DIA-NN only for criteria 1-2)
C++ OpenDIACalibrate (`-panel_targets 40000 -window_panel_targets 300000`) on:
- arm E raw01 and raw02 (the workflow library);
- the PAIRLIB library raw01 and raw02, as a check on runs where the current fit already agreed.

1. **Spread:** spread_s within ±15 % of DIA-NN's robust spread around the written map, on all 4 runs.
2. **Coverage:** DIA-NN 1 % IDs within the written tolerance on every run is ≥ 0.95, and not below the current code's
   value by more than 0.01.
3. **End to end:** self-calibrating OpenDIAWorkflow on both arm-E raws (seed 0).
   - FDP_c at q ≤ 0.01 must be ≤ 0.010.
   - N@FDP_c 1 % must be at least 1718 (raw01) and 1941 (raw02), the values of the current code.

Pass if 1-3 hold. On failure, nothing is retuned and the change is reverted.
