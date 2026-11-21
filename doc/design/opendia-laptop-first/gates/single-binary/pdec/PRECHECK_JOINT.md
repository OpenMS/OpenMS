# PRECHECK_JOINT: would joint modelling (deconvolution, sibling boosting) add information? (21 Nov 2026, diagnostic)

**Question.** Two kinds of joint modelling could lift the missed IDs.

1. **Joint deconvolution.** Fit all co-eluting candidates in an isolation window together, so that a faint precursor whose fragments are shared with a bright neighbour keeps its own residual signal. OpenDIA already has a one-pass version, the interference resolver (`i_explained`, `i_libcos_resid`, `i_dlibcos`, `i_nclaim`). Those four features are among the 45 that PRECHECK_CEIL0 used.
2. **Sibling boosting.** A peptide observed at several charge states should co-elute in all of them. A faint 3+ precursor whose 2+ sibling has a strong peak at the same RT gains support. Decoys have the same charge-state structure as their targets (checked in the r0 PQP: a decoy peptide has decoy precursors at the same charges), so the evidence stays symmetric for target-decoy competition. Entrapment peptides have entrapment siblings, so the entrapment bound polices any correlated error.

RESULT_CEIL0 found that a DIA-NN-labelled ceiling on the current features gains only about +8 %, and that most missed IDs (A_low) have the right window but score like decoys. Joint modelling helps only if it adds information that the current features lack. This precheck measures that on the existing feature tables, with no new runs.

## Setup (fixed before any result)

- **Data:** the PRECHECK_CEIL0 feature tables (`ceil0-20261121/{D,H}_{r0,r7}`, current default, seed 0) with their libraries, and DIA-NN HeLa q ≤ 0.01 IDs from the same level.
- **Script:** `precheck_joint.py`. Its learner and arms are the same as `precheck_ceil0.py`: best-window paired TDC, 3-fold cross-fitting by library fold, HistGBM with 300 rounds and 31 leaves.
- **Sibling features, per window** (label-free; identical rules for targets and decoys; siblings are other precursors of the same PEPTIDE row, meaning the same modified sequence at another charge):
  - `s_n`: siblings in the library;
  - `s_co`: siblings with a window whose apex is within 3 s;
  - `s_best`: the highest tool score among those co-eluting sibling windows, or the run's minimum score if there are none;
  - `s_best_d`: |apex difference| to that window.
- **Arms:**
  - A: the tool score;
  - B: label-free GBM on the 45 features;
  - **Bs:** B plus the sibling features;
  - C: ceiling trained on DIA-NN labels;
  - **Cs:** C plus the sibling features.
- **Interference profile:** computed for the A_low class (the right window is the best one, but q > 0.01), for accepted target best windows, and for decoy best windows. The measures are the median `i_explained` and the share with `i_explained` ≥ 0.5.
- **Sibling profile:** the share of each class with a tool-accepted (q ≤ 0.01) sibling whose best window is within 3 s. For decoys, this means an accepted decoy sibling.

## Decision rule (written before any result)

1. **Sibling boosting is a lever** if all three hold:
   - Bs over B is ≥ +3 % in HeLa IDs at own q ≤ 0.01 at r0, on D and H;
   - Bs over B is ≥ +2 % in N@FDP_c 1 % at r7, on D and H;
   - the r7 lower bound rises by ≤ +0.2 pp.

   If so, opt-in C++ sibling features for the scorer come next, through the usual two-seed check and frozen plan.
2. **Joint deconvolution is a lever** if, on D or H at r0, the A_low share with `i_explained` ≥ 0.5 is ≥ 25 % **and** ≥ 1.5× that of accepted target best windows. The next step is then a trace-level deconvolution precheck: `-traces_out` on a sampled mini library, plus a joint NNLS fit of the co-eluting candidates' library patterns in Python. Otherwise deconvolution is deprioritised, because the missed peaks are not mainly interfered.
3. The ceiling gain Cs over C is reported alongside. If it is < +3 % on both runs, the sibling information is already implicit in what the learner sees.

No code changes follow from this precheck alone.
