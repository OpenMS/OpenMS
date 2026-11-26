# PRECHECK_RTNULL: RT-refinement knockout on used D/H

26 November 2026. Freeze and push this document, scripts and SHA256 manifest before any submission. This is a causal development control, not a new held-out validation and not a declaration of an FDR fix.

## Question and intervention

A source audit found that RT refinement trains on first-pass q <= 0.01 targets, then uses the full fitted ridge model for all library RT centers, including those training targets. Decoys are not training observations. Out-of-fold residuals determine the narrower window but are not the predictions used for training targets' search centers. This might favor training targets. Its FDP effect is unmeasured.

Use the existing `-no_rt_refine` switch as a broad knockout. It removes both the full-fit prediction exposure and the potentially useful second pass/narrower windows, so a change cannot be attributed uniquely to in-sample fitting. A loss of IDs is possible and must be reported, not hidden by relaxing the gate. This intervention uses no species classification or entrapment annotations.

## Fixed cells and inputs

- D and H, original seeds 0,1,2,3,4: ten NORT jobs, 8 threads, 200 GB, 16 h each. No replacement seeds or alternate raw files.
- Baseline: existing `entlib-gbt-20261122/runs/{D,H}_gbt_s{0..4}`. PRECHECK_EVSHARE's representative default gate reported 9/9 files byte-identical across raw03/D/H; this does not mean every baseline seed was rerun in that gate. All baseline files remain read-only. Local baseline precursor tables are fingerprinted in `rtnull/baseline.sha256`.
- Container: `/nfs/research/juan/singj/opendia/evshare-20261126/opendia-evshare-20261126.sif`, SHA256 `c43320a6e0d350efa29d80db2cafe2470bd026edd5a4fbafe6eba6f5792f191a`. Its `--helphelp` was checked for `-no_rt_refine` before freezing. No binary changes.
- Library: `/nfs/research/juan/singj/opendia/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`, SHA256 `76f2a6c22277eaf464a301b446e97b0fd909fbf6b26d79736b1b0337c304f15c`, verified equal to the local evaluation copy.
- Raw paths are fixed in `rtnull/rtnull_hela.sbatch`, identical to the prior D/H runs. Only the extra `-no_rt_refine` flag changes.
- New remote workspace: `/nfs/research/juan/singj/opendia/rtnull-20261126`, verified absent before preparation. Submission refuses existing outputs and repeated submission.
- PXD017703's six held-out runs remain unused. N2 is not part of this control.

## Evaluation and decision, fixed now

Use the unchanged copied `eval_pdec_relabel.py` and existing exact-human step A label audit, fingerprinted with the precheck. Entrapment is accession `instr('_p_target')`, with the existing correction only. No near-homolog exclusions or stratum-specific relabelling. Labels are used only after scoring.

At tool q <= 0.01 report per-seed N, E, lower E/N and original combined E(1+1/r)/N, plus pooled bounds and 90% Poisson intervals exactly as the evaluator supplies. Report median N@FDP_c 1% for each arm and NORT/G0 ratio per set. The library r is unchanged.

A promising calibration candidate requires BOTH sets' original pooled combined bounds <= 1.00% AND BOTH median ID ratios >= 1.00. This is necessary but not sufficient for generalization. Regardless of pass/fail, report the bound changes and ID changes separately. No threshold search, seed dropping, species-dependent rule or tuning after outcomes.

All ten NORT cells must succeed for a pass to be evaluable. Preserve any failure and log. Never substitute a seed. Interrupted infrastructure retries, if ever needed, require a new explicitly documented attempt path and do not erase the first attempt.

## Reproduction and integrity

`rtnull/rtnull_submit.sh` checks the precheck manifest, container and library before submitting. `rtnull/rtnull_eval.sh` evaluates a complete local fetch using fresh symlink layouts without changing baseline tables. The freeze manifest is `PRECHECK_RTNULL.sha256`, relative to this folder. Submission IDs and subsequent outcomes go in separate documents, never this frozen plan.
