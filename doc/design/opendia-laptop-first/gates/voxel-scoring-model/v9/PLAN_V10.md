# Voxel v10: native-sampling fragment evidence with an empirical per-m/z local null (written before any v10 scores were computed)

## Context

- **Base:** D11FWT (PLAN_V6 pass), old population.
- **V9 failure:** PLAN_V9 (D11FNWT) failed on FDP_c, 0.018 / 0.019 (`RESULT_V9.md`). Its uniform-density chance rate was badly mis-calibrated: 0.08 expected matched fragments per away-from-apex decoy scan against 1.38 observed.
- **Pre-plan look** (raw 01 only, descriptive; labels used only to measure; `NATE_SEPARATION_raw01.txt`). Measuring the chance rate empirically per fragment m/z, in the same isolation window, gives these values:
  - decoy windows with apex L ≥ 2 fall to 0.033;
  - the score-matched AUC of `apex_logp` for D against decoys rises to 0.84;
  - entrapment against decoys is 0.55, and all human against decoys is 0.59.
- **Not looked at:** raw 02 features, and any trained score.
- **Risk** stated from the pre-plan look: entrapment against decoys is above 0.5 (0.55 on `apex_logp` and on `cos`). Rule 1.3 is the guard.

## Features E (`native_features_emp.py`; label-free, model-free; same rule for targets and decoys)

Scans, fragments (12-fragment library, calibrated offset, calibrated tolerance) and the 5 columns are identical to PLAN_V9's features N. Only the chance rate changes:

- p_f = (number of scans of the same isolation window, |RT − apex| in (8, 60] s, holding a peak within tolerance of f + 0.5) / (number of such scans + 1);
- p_f is clipped to [1e-3, 0.999].

## Arm

**D11FEWT:**
1. `train_v4.py <raw> D11FE` (11 heuristics + feat12 + E);
2. W rule (`make_w.py`);
3. per-fold tail transform (`calib_tail.py`).

Training, folds, seeds and procedure are the same as for D11FWT.

## Pass rule (both raws; same as PLAN_V7 and PLAN_V9)

1. PLAN_V2 rules 1–4 against BASE. Rule 3 is FDP_c at q ≤ 0.01 ≤ 0.015.
2. N@2 % above D11FWT on both raws.
3. Null-window winner AUC ≤ D11FWT's + 0.005.

**Reported only:**
- N at 0.5 / 1 / 1.5 %;
- exact FDP_c;
- the funnel, including old stage-D misses now accepted;
- the paired bootstrap against D11FWT.

An arm that passes becomes the base.

## Known limits

- **Plan sequence.** This is the second native-evidence plan after a failure. The null was changed after looking at raw 01 feature separation (not at a trained result), and that is stated here.
- **Same limits as PLAN_V9:** co-isolation in 2 Th windows, and possible entrapment asymmetry.
