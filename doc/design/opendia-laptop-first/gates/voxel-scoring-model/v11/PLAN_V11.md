# Voxel v11: native-sampling fragment evidence on unique fragments only (written before any v11 scores were computed)

## Context

- **Base:** D11FWT (PLAN_V6 pass), old population.
- **V10 failure:** PLAN_V10 (D11FEWT) failed on FDP_c, 0.023 / 0.018 (`RESULT_V10.md`). The empirical-null evidence was high for entrapment targets too (median `apex_logp` 5.8–6.0), read (softened) as entrapment peptides matching co-isolated human signal through shared fragments.
- **This plan's idea (user request, last idea before returning to the drawing board):** score only fragments that no other library target in the same isolation window and RT tolerance could explain.
- **Pre-plan look** (raw 01 only, descriptive; labels used only to measure; `NATU_SEPARATION_raw01.txt`):
  - 40.6 % of the 12-fragment slots are shared; median unique fragments is 7 for targets and decoys (mean 6.78 / 7.28);
  - decoy windows with apex L ≥ 2: 0.022 (E: 0.033);
  - score-matched AUC, D against decoys: `apex_logp` 0.851 (E 0.839), `sum_logp` 0.821 (E 0.753), `sig_scans` 0.776 (E 0.701), `persist_logp` 0.780 (E 0.717), `cos` 0.766 (E 0.757), `n_unique` 0.736;
  - score-matched AUC, entrapment against decoys (target ≈ 0.5): `apex_logp` 0.538 (E 0.548), `sum_logp` 0.538 (E 0.491), `cos` 0.545 (E 0.547), `n_unique` 0.528;
  - all human against decoys falls on every column (for example `apex_logp` 0.587 → 0.548).
- **Not looked at:** raw 02 features, and any trained score.
- **Risk** stated from the pre-plan look: the key check (entrapment against decoys ≈ 0.5) is not met. Masking moved it only from 0.548 to 0.538 on `apex_logp` and raised it on `sum_logp`. So the shared-fragment reading of V10 is at most a partial explanation, and the same FDP_c failure is likely. Rule 1.3 is the guard.

## Features U (`native_features_uniq.py`; label-free, model-free; same rule for targets and decoys)

- **Mask:** fragment f (raw library m/z) of precursor p in window w (apex t) is dropped when some other library TARGET q has all of:
  - q is not p and not p's target-decoy partner;
  - q's precursor m/z is inside the isolation window of p;
  - q's predicted run RT (library RT through the run's cal.json knots, `np.interp`) is within the run's RT tolerance (`tolerance_half_width_s`, 84.2 / 83.0 s) of t;
  - q has a fragment within 2 × the calibrated ms2 tolerance of f.
  - Neighbour fragments: the union of the full 6-fragment library (all targets) and the 12-fragment population library (targets).
  - Decoys are not neighbours: they are not real signal and cannot explain a peak.
- **Evidence:** the 5 columns of PLAN_V10's features E (same scans, same empirical per-m/z null, same tolerance), computed on unique fragments only, plus a sixth column `n_unique`.

## Arm

**D11FUWT:**
1. `train_v4.py <raw> D11FU` (11 heuristics + feat12 + U);
2. W rule (`make_w.py`);
3. per-fold tail transform (`calib_tail.py`).

Training, folds, seeds and procedure are the same as for D11FWT. Run by `chain_v11.sh`.

## Pass rule (both raws; same as PLAN_V7, V9 and V10)

1. PLAN_V2 rules 1–4 against BASE. Rule 3 is FDP_c at q ≤ 0.01 ≤ 0.015.
2. N@2 % above D11FWT on both raws.
3. Null-window winner AUC ≤ D11FWT's + 0.005.

**Reported only:**
- N at 0.5 / 1 / 1.5 %;
- exact FDP_c;
- the funnel, including old stage-D misses now accepted;
- the paired bootstrap against D11FWT.

An arm that passes becomes the base. If it fails, this line of work (native-sampling fragment evidence for stage D) stops, as agreed with the user.

## Known limits

- **Plan sequence.** This is the third native-evidence plan after two failures. The mask was set by the user's idea before any look. The pre-plan look shows the key check is not met, and that is stated here.
- **Mask coverage.** Neighbours are limited to the library. Real co-eluting peptides not in the library are not masked. Predicted RT uses library RT and the run's knots, so neighbours with poor RT prediction may be missed or wrongly included.
- **Co-isolation** in 2 Th windows, and possible entrapment asymmetry, as in PLAN_V9 and V10.
