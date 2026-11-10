# Voxel pilot result (raw 01, offline): FAIL, no cell passes

The rule in `PLAN_PILOT.md` was frozen before training (sha256 in `PLAN_PILOT.md.sha256`). The scripts are `make_train_set.py`, `build_voxels.py`, `train_pilot.py` and `evaluate_pilot.py`. The output is `EVALUATION_PILOT.txt`. Training took 13:20 for 24 fits, with a peak RSS of 0.9 GB.

## Result

Held-out folds only. The evaluation sample is the base winning window of 10,068 precursors.

| Cell | AUC present > decoy | AUC entrapment > decoy | N@1 % | N@2 % | ΔN@2 % vs a_T0 [90 %] | Verdict |
|---|---|---|---|---|---|---|
| a_T0 (control) | 0.9989 | 0.4930 | 878 | 1036 | | control |
| b_T0 | 0.9989 | 0.4950 | 869 | 1012 | −24 [−129, +12] | fail |
| c_T0 | 0.9996 | 0.5037 | 476 | 1061 | +25 [−318, +48] | fail |
| d_T0 | 0.9992 | 0.4932 | 780 | 995 | −41 [−155, +7] | fail |
| a_T1 | 0.9983 | 0.4711 | 898 | 1024 | −12 [−64, +28] | fail |
| b_T1 | 0.9974 | 0.4829 | 864 | 1023 | −13 [−69, +16] | fail |
| c_T1 | 0.9988 | 0.4985 | 622 | 1058 | +22 [−229, +45] | fail |
| d_T1 | 0.9983 | 0.4838 | 923 | 1030 | −6 [−78, +37] | fail |
| Base CNN + LDA, same sample (reference) | 1.0000 | 0.5161 | 1028 | 1078 | | |

- No cell has a bootstrap interval for ΔN@2 % above 0. So no cell passes rule 3.
- c_T0 and c_T1 have the highest N@2 % (1061, 1058), but their N@1 % is unstable (476, 622). Their intervals include large losses.
- All T1 cells have an AUC present > decoy below a_T0. So they also fail rule 2.

## Reading

- **Decoy-likeness.** Every voxel cell has an entrapment-over-decoy AUC of 0.47–0.50. The base CNN + LDA has 0.516 on the same sample. So the single voxel model does not learn the decoy artefact that rounds 1–5 found. This is the effect that the design wanted.
- **T1 moves too far.** With decoy-like negatives and synthetic false targets, the entrapment AUC falls to 0.47–0.48 for a, b and d. That is outside 0.5 ± 0.01 on the low side. The model now ranks false targets *below* decoys. This makes the decoy count conservative, but the model may also push down real targets that look like the negatives (not tested).
- **No ID gain.** At FDP_c 1–2 %, no voxel cell reaches the base on this sample (1028 / 1078). The cells use only the voxel input. The base also uses 11 heuristic features (search, boundary, byion), and the byion group alone gave +26 to +36 %.
- **The sample favours the base.** The "present" group is base q ≤ 0.01 winners, and the windows are base winning windows. This makes the base reference optimistic on this sample. It does not change the comparison between the pilot cells.
- **Saturation.** The validation loss goes to 1e-4 or less in most folds. A possible cause is that about 500 anchors are too few to teach the model weak true targets (not tested). Rounds 1–5 showed the same saturation.

## What this does not show

- It does not show that a voxel model cannot beat the base. It shows that these four inputs, with these training sets and 500 anchors, do not beat the voxel control at 2 %.
- The pilot does not include the heuristic features or the b/y ion counts in the model.

## Possible next steps (each needs a new plan before any run)

1. Add the byion counts and the search features to the voxel head, as one change (the base gets most of its gain from them).
2. More positives: semi-supervised rounds on the voxel logit, with the entrapment AUC as a guard.
3. A lower weight for decoy-like negatives (T1 overshoots to 0.47).

The user decides the direction.
