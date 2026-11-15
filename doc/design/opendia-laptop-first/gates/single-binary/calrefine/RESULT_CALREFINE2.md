# PLAN_CALREFINE2 result: fails on raw02; the refined calibration is not adopted

The test ran five seeds per arm with the calibration reused, so only the scorer changes between seeds. The arm-E library
was used, with OpenDIAWorkflow, in `EVAL_CALREFINE2.txt`. The current-code calibration of raw02 reproduces the
self-calibrated workflow run byte for byte.

| | CUR N@FDP_c 1 % (seeds 0-4) | median | REF N@FDP_c 1 % | median | change | FDP_c at q ≤ 0.01, REF |
|---|---|---|---|---|---|---|
| raw01 | 1718, 1687, 1759, 1745, 1709 | 1718 | 1817, 1820, 1813, 1854, 1812 | 1817 | +5.8 % | 0.0057-0.0081 |
| raw02 | 1941, 1999, 1967, 1889, 1908 | 1941 | 1826, 1972, 1880, 1880, 1837 | 1880 | −3.1 % | 0.0034-0.0047 |

- **Criterion 1** (REF ≥ CUR − 1 % on each raw) fails on raw02.
- **Criterion 2** (FDP_c ≤ 0.010 for every REF seed) passes.
- **Verdict:** fail, and nothing changes.
- **CUR FDP_c at q ≤ 0.01:** raw01 0.0048-0.0104 (0.0104 and 0.0095 on two seeds); raw02 0.0023-0.0045.

## Reading
- **The trade-off is consistent over seeds.**
  - raw01 gains about 100 IDs.
  - raw02 loses about 60 IDs at the median, with ranges overlapping.
  - The refined tolerance is narrower on both raws (99 / 100 s against 151 / 124 s). Coverage of DIA-NN IDs falls
    0.967 → 0.959 on raw02, and the run loses true candidates there.
- **Accuracy of the spread estimate is not enough to decide.** The refined spread matches DIA-NN better (ratio 0.96-1.00
  against 1.22-1.39), but 4σ of the true spread leaves about 4 % of true IDs outside the window.
- **The real question is the window width.** The question is what RT window best serves the scorer, not how accurately
  the spread is estimated. Answering it needs a criterion on IDs, not on the spread. That is a design decision outside
  this plan.
