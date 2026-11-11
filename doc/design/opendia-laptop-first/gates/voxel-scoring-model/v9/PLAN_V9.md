# Voxel v9: native-sampling fragment evidence with a spectrum-density chance model (written before any v9 scores were computed)

## Context

- Base: D11FWT (PLAN_V6 pass), old population (`astral-voxel-pilot-20261110/v2`). PLAN_V8 (rtlin population) failed, so the old population stays the reference.
- `v8/stageD-plots/FINDINGS.md` (descriptive; labels used only to define groups):
  - Stage-D misses are real but weak peaks. 75 % hold ≥ 4 of 6 library fragments in the raw apex scan, against 10–13 % of decoy and entrapment winners in the same score range.
  - At **matched** D11FWT score, decoys that reach D's score also hold 4–5 fragments. They get there by chance in dense 2 Th Astral scans.
  - The current features see fragment signal only through 32 RT bins over ± 15 s and through counts that ignore how dense the scan is.
- **Hypothesis.** Scoring each fragment in each raw scan, against the chance of a match given the local peak density, separates real weak peaks from decoys that match by chance. The feature asks how unlikely the observed match count is in this particular scan.

## Features N (`native_features.py`; label-free, model-free; same rule for targets and decoys)

- **Scans.** For each population window, use the MS2 scans of the isolation window that holds the precursor m/z. Take the scan nearest the apex and 3 on each side (revisit 1.65 s), read from `raw0x.mzML` (FileConverter of the same raw).
- **Fragments.**
  - The precursor's 12 library fragments from `library12_pop.pqp`. Decoys use their own decoy m/z.
  - m/z is calibrated with the run's ms2 offset.
  - Tolerance is the run's calibrated half width (`cal.json`, 4.99 ppm on raw 01).
- **Chance model per scan and fragment f.**
  - λ_f = (peaks within ± 5 Th of f, outside the tolerance) × (2 tol_f) / (10 Th − 2 tol_f).
  - p_f = 1 − exp(−λ_f), clipped to [1e-3, 0.999].
  - L = −log10 P(X ≥ k), where X is Poisson-binomial over the fragments' p_f and k is the number of matched fragments. L is capped at 30.
- **Five columns:**
  1. `nat_apex_logp`: max L over the 3 central scans;
  2. `nat_sum_logp`: sum of L over 7 scans;
  3. `nat_sig_scans`: scans with L ≥ 2;
  4. `nat_persist_logp`: −log10 of the chance that as many fragments are matched in ≥ 2 of the 3 central scans;
  5. `nat_cos`: cosine of sqrt(summed central-3 intensities) against sqrt(library intensities) over the 12 fragments.
- **No isolation window.** Windows whose precursor m/z falls outside all isolation windows (119 on raw 01, at the m/z range edges) get zeros.

## Arm

**D11FNWT:**
1. `train_v4.py <raw> D11FN` (11 heuristics + feat12 + N);
2. W rule (`make_w.py`);
3. per-fold tail transform (`calib_tail.py`).

Training, folds, seeds and procedure are the same as for D11FWT. `TRAIN_THREADS=10`.

## Pass rule (both raws; same as PLAN_V7)

1. PLAN_V2 rules 1–4 against BASE. Rule 3 is FDP_c at q ≤ 0.01 ≤ 0.015.
2. N@2 % above D11FWT on both raws.
3. Null-window winner AUC ≤ D11FWT's + 0.005.

**Reported only:**
- N at 0.5 / 1 / 1.5 %;
- exact FDP_c;
- the funnel, including how many old stage-D misses are now accepted;
- the paired bootstrap against D11FWT;
- score-matched AUC of the N features for D against decoys and entrapment against decoys.

An arm that passes becomes the base.

## Known limits (stated before the run)

- **Shared isolation window.** In a 2 Th window, a co-isolated peptide can supply real fragments that are not the precursor's. The density model treats them as chance only through the local peak count.
- **Entrapment asymmetry.** Entrapment targets are real sequences from another proteome. Their fragments may match co-eluting human signal at the same rate as decoys or not. Rule 1.3 (FDP_c) is the guard.
- **Two raws only, and offline Python.** C++ promotion is a separate decision.
