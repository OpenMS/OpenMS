# PLAN_INTERF: fragment-deficit evidence against co-eluting same-composition interference

Written 2026-11-13, before any feature, model or score of this run exists. Root `interf-20261113/`. Option A of
`RESULT_QCAL.md`: keep rule 1c, and add label-free scoring evidence instead.

## Question
FDP_c at q ≤ 0.01 is 0.016 / 0.017 for xic D11FWT, above the 1.5 % limit.
- The excess comes partly from entrapment shuffles that sit on a co-eluting real parent peptide which is missing from
  the library. These match only some of the library fragments.
- Question: do label-free features capturing "some predicted-strong fragments are far below the abundance the others
  imply" bring FDP_c at q ≤ 0.01 to 1.5 % or below on both raws, without losing IDs at fixed entrapment FDP?

## What was looked at before this plan (honest record)
All looks were on raw01, used labels descriptively, and trained no model (`qcal-20261113/interf_look*.py`).
- Rejected candidates:
  - Composition-specificity weights (`interf_features.py`): almost every fragment of a shuffle is unique, and only
    about 6 % of the full ladder is shared with the parent.
  - Extra-fragment shares (`f12_*`).
  - Library-proportion measures (`x_*`).
- Chosen candidate: the intensity-aware deficit.
  - At matched D11FWT score (10 bins above the q ≤ 0.20 cut), d_max has AUC(decoy > human) 0.694 and
    AUC(entrapment > human) 0.702.
  - The raw01 parent-explained entrapments reach 0.740.
- Raw01 therefore informed the feature choice. **Raw02 is the replication and has not been looked at.**

## Arms (same recipe as D11FWT, only the input columns differ)
- `deficit_features.py` computes `featD.npy` per raw. It is label-free, and its columns are defined in the script
  docstring: `d_max`, `d6_max`, `d_lib`.
- **D11FDWT** (primary):
  - `train_v4.py <raw> D11FD`: D11F, that is the 11 heuristics + F12, plus the 3 featD columns. Same seeds, folds,
    training groups, rounds and epochs.
  - Then `make_w.py <raw> D11FD lin_D11FD D11FDW`.
  - Then `calib_tail.py <raw> D11FDW D11FDWT`.
- **D11F control:** `train_v4.py <raw> D11F` in this root, then W and T into D11FWT_rerun. This checks that
  training reproduces the xic D11FWT within run-to-run noise. It is reported, not ruled on.
- Inputs are symlinked from `xicids-20261112/<raw>/`: population, heuristics, feat12, 12-fragment arrays, train
  windows, and BASE byion features. Nothing is rebuilt.

## Pass rule, D11FDWT on both raws (PLAN_XIC_IDS rules, unchanged)
1. Against the rebuilt BASE (`byion/features.tsv`, as in PLAN_XIC_IDS):
   - The bootstrap 90 % interval of N@2 % − BASE lies above 0.
   - The interval of N@1 % − BASE is not entirely below 0.
   - **FDP_c at q ≤ 0.01 is ≤ 0.015.**
   - AUC(entrapment > decoy winner) is ≤ the BASE AUC.
2. N@FDP_c 2 % is > 1716 (raw01) and > 1736 (raw02).
3. The null-window winner AUC is ≤ 0.4983 (raw01) and ≤ 0.4986 (raw02).

Recorded but not ruled on:
- N@FDP_c 0.5/1/1.5/2 % against xic D11FWT;
- E, D and rho at q ≤ 0.01, 0.05 and 0.10 (`qcal-20261113/diagnose.py`);
- the parent-explained share (`parent.py`).

## Run
`chain.sh`, run locally: features, then 2 arms × 2 raws of training (about 2 min each), then W, T, evaluation and
null AUC, for about 20 min in total. One run, with no retuning after a failure. If D11FDWT fails, the result is
reported as it is, and any change goes into a new plan.
