# PLAN_V19: native-sampling peak geometry replaces the search-run summaries (D11FGWT)

Written after the raw01 feature look and before raw02 features or any V19 model/score. Base: D11FWT.

## Finding that motivates it
The D11F inputs `transition_support`, `log_matched_scans` and `log_peak_width` come from OpenDIASearch's streaming
signal run (`updateCandidate_`): a run ends only after a 30 s gap, so width and scan count measure the run, not the
peak (median width 46-66 s against DIA-NN FWHM 3-4 s; Spearman with DIA-NN FWHM 0.02), and support is the maximum
fragment count in any scan of the run. Model attribution on DIA-NN(our-library)-only IDs (`gates/diann-ourlib`) put the
largest nonlinear penalties on these three inputs.

## Features G (`build_peak.py`, `peak_core.py`; label-free, model-free, same rule for targets and decoys)
Per window: the MS2 scans of the precursor's isolation window, scan nearest the search apex and 4 on each side; the 6
search-library fragments (calibrated offset and tolerance, as V9). Trace = per-scan summed matched intensity.
Native apex = trace maximum within +-2 scans of the search apex.
- col 0 -> fragments matched in >= 2 of the 3 scans centred on the native apex;
- col 1 -> log1p(contiguous scans above half maximum around the native apex);
- col 2 -> log1p(FWHM in s, linear half-maximum crossings, trace edge if never below half).
All other inputs, rows and window identities are unchanged; 17 inputs as before.

Raw01 look (labels used only to summarise): native FWHM vs DIA-NN FWHM Spearman 0.56 (old 0.02). AUC vs decoy winners:
native support DN 0.84 / ENT 0.54 (old 0.91 / 0.55); native points 0.61 / 0.52 (old 0.44 / 0.52); native FWHM
0.66 / 0.53 (old 0.41 / 0.48). Risk stated: native support alone is weaker than the old support.

## Arm
D11FG = D11F with columns 0-2 replaced, trained by `run_v19.py` (imports `train_v4.py` unchanged; same seeds, rounds,
thresholds, validation hash, objective, 10 threads; in-runner D11F control must reproduce within 1e-3 or stop).
Then `make_w.py` (D11FGW) and `calib_tail.py` (D11FGWT), unchanged. One run per raw, no tuning.

## Pass rule (both raws; unchanged from PLAN_V7 / V9-V11 / V18)
1. PLAN_V2 rules 1-4 against BASE (N@2 % interval above 0; N@1 % interval not entirely below 0; FDP_c at q <= 0.01
   <= 0.015; entrapment-over-decoy winner AUC <= 0.5313 / 0.5348).
2. N@2 % above D11FWT on both raws (1716 / 1736).
3. Null-window winner AUC <= D11FWT + 0.005 (0.4933 / 0.4936).

## Reported only
N at 0.5 / 1 / 1.5 %; paired bootstrap vs D11FWT; q <= 0.01 human / entrapment / paired-decoy counts, gained and lost
precursors and peptides vs D11FWT (`gains.py`); old stage-D gains and E losses; overlap with DIA-NN on our library
(raw01). Absolute entrapment and decoy increases are stated as costs even if the rule passes.

## Known limits
- The candidate set and search apex are unchanged; this only re-measures each candidate window.
- The 6-fragment trace sum is dominated by the strongest fragment.
- Raw02 is not a pristine holdout.
