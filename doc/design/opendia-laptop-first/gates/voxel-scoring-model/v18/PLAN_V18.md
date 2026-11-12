# PLAN_V18: charge-sibling co-elution feature (D11FSWT)

Written before any V18 feature, model or score was computed. Base: D11FWT (old population, raws 01/02).

## Motivation
V18 look (`../v18_look/`, design docs `v18/LOOK_V18.md`): co-elution of another charge state of the same peptide is
entrapment-neutral (ENT vs DEC AUC 0.508 / 0.504) and separates D from decoys (0.632 / 0.656). It is measured in
another precursor's data, unlike all previous evidence arms.

## Features S (`build_sib.py`, `sib_core.py`; label-free, same rule for targets and decoys)
For every candidate window, over other precursors of the same library peptide (decoys: the same decoy peptide) in the
population, with z = per-fold decoy-winner median/MAD normalised existing held-out D11F full score (winners by lin_D11F,
the W rule):
1. has_sib: the peptide has another precursor in the population;
2. has_near: a sibling window has apex within 3.0 s of this window's apex;
3. sib_z_near: max sibling z within 3.0 s, clipped to [-5, 20]; 0 when has_near = 0.
Siblings share a family and so a fold: held-out rows only use held-out D11F scores. Not used: entrapment, DIA-NN.
Threshold 3 s and no others were fixed from the raw01 look. No tuning.

## Arm
D11FS = D11F inputs (11 heuristics + feat12) + S, trained by `run_v18.py`, which imports `train_v4.py` unchanged and
repeats its D11F loop (same seeds, rounds, thresholds, validation hash, objective, 10 Torch threads). The same runner
first retrains the D11F control and stops unless it reproduces scores_D11F / lin_D11F within 1e-3. Then
`make_w.py` (D11FSW) and `calib_tail.py` (D11FSWT), unchanged. One run per raw, no seed or threshold search.

## Pass rule (both raws; unchanged from PLAN_V7 / V9 / V10 / V11)
1. PLAN_V2 rules 1-4 against BASE: 90 % interval of N@2 % minus BASE above 0; 90 % interval of N@1 % minus BASE not
   entirely below 0; FDP_c at q <= 0.01 <= 0.015; entrapment-over-decoy winner AUC <= BASE (0.5313 / 0.5348).
2. N@2 % above D11FWT on both raws (1716 / 1736).
3. Null-window winner AUC <= D11FWT's + 0.005 (0.4933 / 0.4936).

## Reported only
N at 0.5 / 1 / 1.5 %; paired bootstrap vs D11FWT; q <= 0.01 human / entrapment / paired-decoy counts with gained and
lost precursors and human peptides vs D11FWT (`gains.py`); old stage-D gains and stage-E losses; funnel.
The user's aim is more true targets without more entrapments or decoys: absolute entrapment and decoy counts are
reported next to the pass rule, and an increase is stated as a cost even if the rule passes.

## Known limits
- Stacking: S on training rows comes from D11F models that saw the test fold as training data (indirect, symmetric
  for targets and decoys). Test-fold rows use only held-out scores.
- Correlated errors: a false charge state can support another false charge state of the same peptide. Decoys have the
  same mechanism (89 % of decoy winners have a library sibling), so TDC and entrapment measure it.
- Gains may be extra charge states of peptides already found; peptide-level counts are reported.
- Raw02 is not a pristine holdout.
