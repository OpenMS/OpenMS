# PLAN_RTWIN2: RT window width chosen from the run itself (cross-fitted)

Written 16 November 2026, before the implementation and before any run. It replaces the label-free density rules,
which were undefined (`PRECHECK_CROSSING.md`).

**User requirement.** The width must come from the file being analysed. A single file has to work, and so do
heterogeneous files. No width may be tuned on other datasets.

## Design
1. **Search wide.** OpenDIAWorkflow searches inside the calibrated RT tolerance T. The current rule is 4 × the spread;
   it is wide on raw 01. There is no change up to scoring.
2. **Score.** The scorer runs as now: features, three folds, network, W rule and tail calibration.
3. **Choose the width per fold, cross-fitted.**
   - For fold k, and for each candidate half width w on a grid of T × {0.40, 0.45, …, 1.00} (13 values), restrict every
     precursor to the windows whose |apex − calibrated library RT| ≤ w.
   - Take each precursor's kept window as its W-rule window (largest linear part, as now) among its windows inside w.
     Then compute the number of target
     precursors at pooled paired TDC q ≤ 0.01 on the **other two folds** only, with their decoys.
   - The grid value with the largest count is w_k. Ties go to the larger w.
   - Fold k is then reported with w_k, so a fold's width never sees its own targets or decoys.
4. **Scores and q-values.** The scores are not retrained. A window outside w_k is dropped for its precursor; if a
   precursor has no window left inside w_k, it is dropped too.
   - Tail calibration and pooled paired TDC are recomputed on what remains.
   - Each fold's decoys are restricted with that fold's own width, the same as its targets.

This is the same competition the run reports. It uses only target/decoy status, never entrapment labels. The width is
chosen per run (per fold) and costs no extra raw pass.

## Arms (OpenDIAWorkflow, arm-E library + full 12-fragment scoring library `library12_full.pqp`, seeds 0-4)
- **CUR:** self-calibrating, with no width selection. This is the current code with the full scoring library.
- **RTWIN:** self-calibrating, with the width selection above.

Both arms come from one run per seed, because RTWIN is a post-scoring step. With `-rt_width_selection` the tool writes
the selected table as `<run>.precursors.tsv` and the unselected one as `<run>.precursors.full_window.tsv`.

**Runs.**
- Raws 01 and 02 are the design runs.
- **Raw 03 is held out.** It is evaluated once, after the verdict on raws 01/02 is written. Raw 03 as the held-out run
  was part of the plan the user said to go ahead with on 16 November 2026.

## Criteria (entrapment labels only in evaluation; eval_run.py)
On raws 01 and 02 separately, both of the following must hold:
1. **Median N@FDP_c 1 %:** RTWIN ≥ CUR − 1 %.
2. **FDP_c at q ≤ 0.01:** ≤ 0.010 for every RTWIN seed.

The pass also needs RTWIN ≥ CUR + 2 % on at least one raw: the step must help somewhere, not only be neutral.
Raw 03 is then reported against the same two criteria (descriptive, held out).

Pass: RTWIN becomes the default in OpenDIAWorkflow. Fail: it stays an option, default off. Nothing is retuned.
