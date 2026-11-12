# XIC peak picking (no 30 s gap), raw01 look

12 November 2026. Code: OpenMS-r1-sparse-experiment 11522c7a6c (`OpenDIA -peak_picking xic`, default unchanged).
Plan frozen before any output: `PLAN_XIC_LOOK.md` (sha256 in `PLAN_XIC_LOOK.md.sha256`). Read-outs: `LOOK_raw01.txt`,
extra splits `LOOK_raw01_extra.txt`. No scoring change yet.

## What changed
The search still uses the same streamed grid matching (same tolerances, calibration and isolation lanes). Previously a candidate
was a signal run closed after a 30 s silent gap, with apex from summed intensity. Now every fragment's matched intensity is
kept at every native scan inside the calibrated RT window, and when the stream passes the window the precursor's traces are
peak-picked from their shape: local maxima merged unless the valley is below half the lower one, boundaries at the valley or
5 % of apex, FWHM, and fragment co-elution (sum of per-fragment correlations with the other fragments). Top 5 peaks per
precursor by co-elution. No gap or width constant is left; the only window is the calibrated RT window.

## Read-outs (go rule: all four pass)

| | legacy (30 s gap) | xic | rule |
|---|---:|---:|---|
| Wall, raw01 search | 2:46 (2:47–3:20 earlier) | 3:17 | ≤ 1.5x: **pass** |
| Max RSS | 4.34 GB | 4.45 GB | ≤ 12 GB: **pass** |
| Legacy output byte-identical | | IDENTICAL (features and peak table) | **pass** |
| DIA-NN IDs with a candidate within 3 s | 90.5 % | **92.6 %** | xic ≥ legacy: **pass** |
| Within 10 s | 96.4 % | 96.7 % | |
| Width vs DIA-NN FWHM, Spearman | 0.025 | **0.557** | ≥ 0.5: **pass** |
| Median width on DIA-NN IDs | 46–54 s run length | FWHM 3.06 s (DIA-NN 3.50 s) | |

- DIA-NN-only IDs (DN, 643): a candidate within 3 s rises from 69.7 % to 77.1 %. For IDs both accept (1692) it stays at 98.4 %.
- The best-co-elution peak (rank 0) is the one nearest DIA-NN's RT for 94 % of matched IDs. For DN the rank-0 peak is within
  3 s for only 61.7 % (BOTH 97 %). So the right peak is usually there, but ranking DN peaks needs the scorer.
- Rows: 2.68M (legacy 2.60M), 579K precursors with a peak (legacy 515K), about 10 maxima per RT window before the top-5 cut.
- Untrained separation (AUC vs decoy rank 0): co-elution DN 0.869, BOTH 0.992, entrapment 0.522. That is single-feature only;
  legacy best support gets 0.919 on DN with a 30 s run, so co-elution alone is not yet the better single score.
- A single-score TDC on rank-0 co-elution alone gives 258 at q ≤ 0.01 with 0 entrapments. This is only a sanity check.
- FWHM floor: 64–68 % of rank-0 peaks have one scan above half maximum (FWHM 1.65 s = one revisit). Weak peaks on a
  ~1.65 s revisit are at the sampling limit. An interpolated or fitted FWHM could help later.

## Next
Go: rebuild the population and stage-2 features on xic candidates, then retrain D11F on it under an unchanged pass rule.
This needs a new population (the old one is keyed to legacy candidate rows), raw02 search, and a frozen PLAN before scoring.
