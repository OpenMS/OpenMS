# PLAN_GAP: three label-free evidence blocks for the right-peak-but-rejected IDs

Written 2026-11-14, before any score with these features was computed. Root:
- local: `OpenMS-r1-benchmark-20261107/gap-20261114/`
- codon: `/nfs/research/juan/singj/opendia/gap-20261114/home`, bound onto the PAIRLIB path.

## Question
`RESULT_GAP.md` found that 626 / 644 DIA-NN 1 % IDs (raw01 / raw02) are lost at the scoring step:
- OpenDIA's winner is DIA-NN's peak, but the winner ranks below the N@FDP_c 1 % cut.
- These IDs are weak (fewer fragments with signal), and they sit far below the cut.

User request: try three fixes, each as one block of label-free features added to the arm E scorer (D11FDWT):

| Arm | Block | Script | Columns |
|---|---|---|---|
| D11FDA | A: fragment mass consistency and dominance | `features_A.py` | a_ppm_llr, a_ppm_sd, a_top_share, a_top_excess |
| D11FDG | G: run-specific fragment intensity correction learned from the calibration anchors, then 12-fragment evidence recomputed with it | `features_G.py` | g_libcos, g_support6, g_coel6, g_dmax |
| D11FDI | I: explaining away shared fragment signal; claimant precursors' expected contribution is removed, and evidence is recomputed on the residual (step 2 of `DESIGN_INTERFERENCE_GRAPH.md`, a one-pass fit) | `features_I.py` | i_explained, i_libcos_resid, i_dlibcos, i_nclaim |
| D11FDAGI | all three blocks | | 12 columns |

Definitions are in each script's header.
- Every block is computed for targets and decoys by the same rule, with no labels.
- Every constant either comes from the run (cal.json m/z model and tolerance, peak extent, noise floor) or is inherited from the existing feature scripts: central bins 12:20, extraction window ±17.5 ppm, co-elution threshold 0.5.

## Recipe
- **Data:** arm E's run root, unchanged: library, 12-fragment PQP, population, voxels, xic12, heuristics, feat12, featD and calibration. Files are hard-linked read-only, and scores and lin files are not linked.
- **Training:** `train_gap.py` is `train_v4.py` plus the three optional blocks, appended after featD in the order A, G, I. Diff in `scripts/`.
- **Control:** D11FD is retrained in the same job. It must reproduce arm E's D11FDWT, which is checked in `CONTROL.txt`.
- **Post-processing, unchanged:** make_w (W), calib_tail (T), evaluate_v2, null_auc, diagnose and pooled_rho.
- **New:** `compare_gap.py` gives the paired bootstrap against the control D11FDWT, with evaluate_v2's resampling and seed.
- **Compute:** one codon job (`run_gap.sbatch`, 12 CPUs, 96 GB).

## Pre-plan look (disclosed)
- Before this plan, I computed the three blocks on **raw01 only**.
- I compared each column descriptively, using evaluation labels, between two groups: the 626 rejected IDs against score-matched decoy winners, and 128 entrapment winners in the same score band against their own matched decoys (`check_featA.py`).
- Best column per block, AUC rejected>decoy / entrapment>decoy:
  - A a_ppm_llr: 0.750 / 0.551
  - G g_libcos: 0.730 / 0.539
  - I i_libcos_resid: 0.769 / 0.507
- G's anchor model reduces CV error by 6.2 % on raw01.
- No scorer was trained, and no definition was changed after the look except one, recorded here: G's anchor-match tolerance was first a fixed 2 × 2.66 s and was replaced by the data-derived value (peak extent × median half width = 1.81 s). The look was rerun, and the AUCs were identical.
- Raw02 has not been looked at with these features.

## Pass rule (per arm, D11F…WT, both raws, against the control D11FDWT of this run)
1. **Specificity:**
   - FDP_c at q ≤ 0.01 ≤ 0.015 on both raws.
   - AUC(entrapment > decoy winner) ≤ control + 0.005.
   - Null-window winner AUC ≤ control + 0.005.
2. **Gain:**
   - The ΔN@FDP_c 2 % 90 % interval lies above 0 on both raws.
   - ΔN@FDP_c 1 % has a point estimate > 0 and an interval not entirely below 0 on both raws.
3. **Calibration:**
   - Pooled rho at q ≤ 0.10 ≤ 1.3.
   - The pooled rho 90 % interval at q ≤ 0.05 includes 1.

**Choice:** among the passing arms, the one with the highest summed N@FDP_c 1 %. If AGI passes and beats the best single block by less than its own ΔN@1 % interval half-width, the simplest passing single block is named as well.

**Expectation, written now:** each single block +50 to +150 at N@FDP_c 1 % per raw; I ≥ A ≥ G. Not part of the rule.

## Recorded, not ruled on
- N@FDP_c 0.5 / 1 / 1.5 / 2 %.
- The q ≤ 0.01 N.
- rho per q cut.
- Funnel recovery: how many of the 626 / 644 right-peak-rejected IDs are accepted, and how many previously accepted are lost (`funnel.py` with the arm's score, after the run).
- Wall time and memory per step.

## Caveats written in advance
- One run per arm and no retuning after a failure.
- Raw02 is replication, but the raw01 look above means raw01 is not a clean test for these blocks.
- I is a one-pass approximation: a claimant's abundance comes from its own unshared fragments, without the iterative joint fit. A pass would justify the full graph model. A fail does not rule it out.

## Amendment 1 (2026-11-14, after the chain stopped at `compare`; scores and EVALUATION.txt existed and were seen)
- **Cause:** `compare_gap.py` exec'd evaluate_v2.py's header, which redefines `ARMS` and so overwrote the arm list. It
  then looked for the evaluate_v2 default arms (`scores_V0.npy`) and stopped. No score, feature or rule is involved.
- **Fix:** the arm list is kept in `GAP_ARMS` and restored after the exec. `chain_gap_eval.sh` reruns only the
  remaining evaluation steps of chain_gap.sh (compare, null_auc, diagnose, pooled_rho, control check), unchanged.
- **Disclosure:** EVALUATION.txt (N@FDP_c, q <= 0.01 FDP_c, entrapment AUC per arm) was read before this amendment.
  The pass rule, arms and choice are unchanged.
