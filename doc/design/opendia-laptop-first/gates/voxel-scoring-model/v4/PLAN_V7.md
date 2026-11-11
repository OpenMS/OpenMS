# Voxel v7: cross-run evidence (written before any v7 scores were computed)

**Context:**
- Base: D11FWT (PLAN_V6 pass).
- `DIAG_CROSSRUN.txt` (diagnostic on existing D11FWT scores; DIA-NN and entrapment labels used only to measure):
  - The other run's D11FWT score, at a window within 10 s after RT alignment, separates the right windows of stage-D misses from the top-1 % decoy winners. AUC 0.736 / 0.749 (raw 01 / 02), against 0.537 / 0.633 for the run's own score.
  - 78–80 % of stage-D misses have such a window. So do 46–48 % of the top-1 % decoys and 20 % of all decoy winners.
- The model must not train on another model's scores. So v7 uses the other run's **raw features**, not its score.

## Features X (`xrun_features.py`; label-free, model-free; same rule for targets and decoys)

For each population window w of precursor p in run a, take the other run's window of the same precursor whose apex, mapped to run a, is nearest to w's apex and within 10 s.

- **RT map b → a:** 20 equal-count binned-median knots of `run_rt` of the calibration anchors (`anchors.tsv`, used = 1) shared by both runs (441 anchors; residual median 0.39 s).
- **19 columns:**
  - that window's 11 heuristics and 6 feat12 values;
  - |mapped apex − w apex| in s;
  - a presence flag.
- **Absent window:** presence 0, |drt| 10 s, and the 17 feature columns set to the median of the present rows.
- **Check:** presence is 0.570 / 0.566 for target / decoy windows on raw 01, and 0.569 / 0.564 on raw 02.

## Arm

| Arm | Construction |
|---|---|
| D11FXWT | `train_v4.py <raw> D11FX` (11 heuristics + feat12 + X), then the W rule (`make_w.py`), then the per-fold tail transform (`calib_tail.py`). Same training, folds, seeds and procedure as D11FWT. |

## Pass rule

Both raws; same as PLAN_V6, with the control now D11FWT.

1. PLAN_V2 rules 1–4 against BASE.
2. N@2 % above D11FWT on both raws.
3. Null-window winner AUC ≤ D11FWT's + 0.005.

**Reported only:**
- N at 0.5 / 1 / 1.5 %;
- exact FDP_c;
- the funnel;
- the paired bootstrap against D11FWT.

An arm that passes becomes the base.

## Known limits (stated before the run)

- **Scoring couples runs.** A run's score depends on the other run. In C++ this would mean a multi-run scoring step, a change of architecture. That is the user's decision before any promotion.
- **Only one other run.** Two raws give one other run, so these results do not show how X behaves with more runs, or with a single run.
