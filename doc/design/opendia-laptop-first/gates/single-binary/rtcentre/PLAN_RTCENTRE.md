# PLAN_RTCENTRE (frozen before any run with the change)

## Change
In the one-pass calibration (`CalibrationCollector::calibrate`, no window), stage 3 is centred on the **stage-1 map
(map1)** instead of the stage-2 EM map. Library RTs outside the stage-1 anchor support
`[stage1.front().first, stage1.back().first]` are searched to the run edge on that side, with no RT term, exactly as
in the windowed pass (`finish_(..., support_low, support_high)`).
- The stage-3 half-window (max(4σ1, 30 s)) and every other step are unchanged. em_map is still computed and reported.
- The windowed pass is not touched. Astral (raws 01-03 all calibrate in 2 passes) must be byte-identical.
- No new constants and no option: this is the default path.

Motivation: PRECHECK_RTCENTRE.md. On HeLa the EM map is 16-29 s late at the gradient start, while map1 is unbiased.

## Runs
- Binary: the OpenDIAWorkflow commit with the change, built in `OpenMS-r1-full-ninja-build`. Baseline: 0e9e4fc3a5
  (copy at `~/.jcode/scratch/dia2/bin/OpenDIAWorkflow.base`).
- Design run D: HeLa 200 ng `20200505_..._2737.d`, library `200623_hela_evosep_fractions_library_ptypic_decoy.pqp`.
  Seeds 0-4, NEW and BASE, 8 threads, local.
- Held-out run H: HeLa 200 ng `20201207_tims03_..._high_speed_21min_..._22474.d`, same library. Seeds 0-4, NEW and
  BASE. Evaluated only after the verdict on D and A is committed.
- Astral A: raw01 seed 0, NEW. Must be byte-identical to BASE (precursors.tsv, run precursors.tsv, calibration.json).

## Metrics (no entrapment on diaPASEF; DIA-NN 2.2 on the same library is evaluation-only)
- N: target precursors at the tool's own q ≤ 0.01 (`precursors.tsv`).
- Agreement: share of N in DIA-NN's q ≤ 0.01 list (a specificity proxy).
- Map: |calibrated library RT − DIA-NN RT| p95 for DIA-NN's IDs.

## Pass rules (all must hold)
1. D: median N(NEW) ≥ median N(BASE) × 1.02, and min N(NEW) > max N(BASE).
2. D: median agreement(NEW) ≥ median agreement(BASE) − 0.005.
3. D: map p95(NEW) < map p95(BASE).
4. A: byte-identical.
5. H (after the verdict): median N(NEW) ≥ median N(BASE) × 1.00, and median agreement(NEW) ≥ median agreement(BASE) − 0.005.
   If H fails, the change is reverted even if D passed.

On failure: revert, no retuning.
