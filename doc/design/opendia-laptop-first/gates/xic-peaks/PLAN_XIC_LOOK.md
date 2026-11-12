# XIC peak picking, raw01 look (descriptive, no scoring change)

Written 2026-11-12 before any xic output exists. Frozen by sha256 (PLAN_XIC_LOOK.md.sha256).

## What changed
`OpenDIA -peak_picking xic` (new; default `run` unchanged). Same streamed grid matching (same tolerances,
same calibration, same isolation lanes). Instead of closing a signal run after a 30 s gap, every fragment's matched
intensity is kept at every native scan of the precursor's isolation window inside its calibrated RT window
(+-84 s). Once the stream passes the window the precursor's native traces (zeros where nothing matched) are
peak-picked (`pickXicPeaks`): consensus = sum of fragments, [1 2 1] smoothing, local maxima merged while the valley
is above 0.5 of the lower maximum, boundaries at the valleys or 5 % of apex. Apex, area, FWHM from the raw
consensus; co-elution = sum over fragments of max(0, r) with the sum of the other fragments; ranked by co-elution,
top 5 per precursor. No gap, width or RT constant other than the calibrated RT window. `-min_transitions 1`.

## Runs (raw01, same .raw, library, cal.json as rtlin-20261111/raw01)
1. xic: `-peak_picking xic -xic_max_peaks 5 -min_transitions 1`, `/usr/bin/time -v`.
2. legacy check: the new binary with the frozen rtlin search options; `cmp` against rtlin-20261111/raw01/pipeline/search.features.tsv.

## Read-outs (labels used only here; DIA-NN = diann-ourlib-20261112, q <= 0.01 targets, RT/FWHM in min -> s)
- R1 speed/memory: wall and max RSS of xic vs legacy (3:20, 4.32 GB).
- R2 legacy byte-identical (must be IDENTICAL, else the change is not safe to keep).
- R3 apex agreement: share of DIA-NN IDs with any candidate apex within 3 s / 10 s of DIA-NN RT, xic vs legacy; and
  for xic the rank-0 candidate.
- R4 FWHM: Spearman of xic FWHM vs DIA-NN FWHM on the candidate nearest DIA-NN RT within 3 s (legacy width 0.02,
  V19 native 0.56).
- R5 untrained separation: AUC of the rank-0 co-elution, DIA-NN-only (DN), both (BOTH, D11FWT q<=0.01 too),
  entrapment best (ENT) each vs decoy rank-0 (DEC). Legacy reference: best transition_support per precursor.
- R6 untrained single-score TDC on rank-0 co-elution (target-decoy pairs, q<=0.01): N, entrapment count, FDP_c.
  Descriptive only: it shows what the new peaks give before any model.
- R7 rows and peaks per precursor.

## Go rule for the next step (rebuild population + features on xic candidates)
Go if all: R2 IDENTICAL; R1 wall <= 1.5x legacy and RSS <= 12 GB; R3 xic any-candidate within 3 s >= legacy within 3 s;
R4 Spearman >= 0.5. Otherwise report and fix the picker before going further. No changing these thresholds after looking.
