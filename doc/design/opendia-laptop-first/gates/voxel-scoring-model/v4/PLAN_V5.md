# Voxel v5: fragment mass accuracy (P) and consensus co-elution (C) for the single decoy-whitened network (written before any v5 training)

**Context:**
- Base is D11FW (PLAN_V4B pass, `RESULT_V4P2.md`). PLAN_V4B names the arm that passes as the new base.
- `SPARSE_EVIDENCE.txt` / `SPARSE_GRID.txt` (diagnostics, DIA-NN reference only):
  - The stage-D misses are not "few strong fragments + strong MS1". Only 3 / 4 precursors have that pattern.
  - Below the cut, every tier of strong fragments has a target-to-decoy ratio of 1.1–1.6.
- `MASS_ACCURACY.txt`: the number of fragments within 3 ppm separates stage-D right windows from the top 1 % of decoy winners (AUC 0.74 / 0.81). No current feature carries this (largest correlation with the 11 features: 0.40).
- `ROBUST_COELUTION.txt`: median-consensus co-elution separates stage D from below-cut decoys (AUC 0.71 / 0.73), but not from the top decoys (0.46 / 0.49).

## Features

`features_v5.py`, label-free, from the existing arrays:
- **P:**
  - `p_n_tight`: of the 6 original fragments, those with signal at the apex (±4.7 s) and |mean ppm − calibrated offset| ≤ 3;
  - `p_ppm_abs_med`;
  - `p_tight_share`.
- **C** (12 fragments, each scaled to its own maximum, so weak fragments count as much as strong ones):
  - `c_med_top3`: mean of the 3 highest correlations with the per-bin median consensus;
  - `c_med_k`: fragments with that correlation ≥ 0.8.

## Arms

All arms use the D11FW settings: `train_v4.py` with α = 0, same seeds, population, folds and rounds, and the PLAN_V4B W window rule on the final scores (`make_w.py`).

| Arm | z |
|---|---|
| D11FPW | 11 + F (6) + P (3) |
| D11FPCW | 11 + F (6) + P (3) + C (2) |

## Pass rule (each arm, both raws)

1. PLAN_V2 rules 1–4 against BASE:
   - the ΔN@2 % interval is above 0;
   - the ΔN@1 % interval is not entirely below 0;
   - FDP_c at q ≤ 0.01 is at most 0.015;
   - the entrapment AUC is at most BASE's.
2. N@2 % is above D11FW on both raws (paired bootstrap against D11FW reported).
3. Null-winner AUC ≤ D11FW's + 0.005.

**Reported only:**
- the funnel against DIA-NN at 1 % and 2 %;
- N at 0.5 / 1 / 1.5 %;
- exact FDP_c;
- the paired bootstrap against D11.

An arm that passes becomes the base.

## Not in this plan (later, each with its own plan)

- Per-fold tail calibration (1 % instability on raw 02).
- A lower population threshold for sparse precursors (2 tight, co-eluting fragments).
- RT recalibration for stage B.
