# PRECHECK_IMPICK: would RT/IM (2D) peak picking pick the right peak more often? (19 Nov 2026, diagnostic)

**Data.**
- r ≈ 1.56 HeLa D seed 0, current default (`-features_out`, `gap7/ev/D_new_s0`).
- Truth, for evaluation only: DIA-NN 2.2 RT and 1/K0 for its 25,881 non-entrapment q ≤ 0.01 IDs on the same raw file.
- A window is "right" if its apex is within 6 s of DIA-NN's RT.
- Scripts: `precheck_impick.py` and `precheck_impick2.py`.

## Where the DIA-NN IDs end up

| outcome | Q1 low | Q2 | Q3 | Q4 high | all |
|---|---|---|---|---|---|
| accepted, right peak | 0.197 | 0.328 | 0.532 | 0.836 | **0.473** |
| missed: right peak kept but scored low | 0.453 | 0.380 | 0.280 | 0.105 | **0.305** |
| missed: right peak never a candidate | 0.222 | 0.221 | 0.147 | 0.042 | **0.158** |
| missed: right peak a candidate, another kept | 0.126 | 0.068 | 0.039 | 0.014 | **0.062** |
| accepted, other peak | | | | | 0.003 |

This corrects an earlier summary. Most wrong-peak losses (16 % of IDs) are peaks the search never produced. Only 6 % of IDs are lost by choosing the wrong one among existing candidates.

## Choosing among candidates: mobility does not separate right from wrong

Across the 1,644 choice errors, the kept window was compared with the best right candidate:

| feature | right is better in |
|---|---|
| mobility centroid offset | 0.53 of pairs |
| central mobility share | 0.49 |
| mobility co-elution | 0.48 |
| \|library RT error\| | **0.82** |
| MS2 area | 0.63 |

- Only 11 % of the kept wrong windows have a mobility centroid more than 2 bins (0.015 1/K0) off centre, against 8 % of right windows.
- The tool already keeps a right window for 92.4 % of IDs that have one. Picking by mobility alone would manage 44-54 %.

## The right peak that is never a candidate (15.9 % of IDs)

| cause | share |
|---|---|
| DIA-NN RT outside the calibrated RT window (±35 s) | 1.0 % |
| DIA-NN 1/K0 outside the calibrated 1/K0 window (±0.052) | 0.1 % |
| inside both windows, but the 5-candidate cap was hit | 6.7 % |
| **inside both windows and below the cap** | **92.1 %** |

- The 1/K0 window is about 6× wider than the real error: |1/K0 error| is p50 0.009 and p95 0.031, against a ±0.052 window.
- So the search's 1-D XIC pools every ion inside ±0.052. A co-eluting precursor at a different mobility inside that
  window can mask or distort the true peak.
- This is the one place where 2-D (RT × IM) picking, or a narrower mobility window at search, could add candidates.
  Its size is not yet measured. The other possibility is that these peaks are simply too faint: they sit mostly in Q1-Q2.

## Conclusion
- **Choosing between peaks:** 2-D picking is not supported as a fix. Mobility does not separate right from wrong
  candidates, and RT already decides correctly most of the time.
- **Missing candidates (16 % of IDs):** 2-D picking is plausible and should be tested next. A cheap first test is a
  search with a narrower 1/K0 window, with candidate recall at DIA-NN RT as the label-free-safe measure.
- **The largest loss is unchanged:** right peak, low score, 30.5 % of IDs, mostly at low abundance.
