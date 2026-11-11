# Voxel v4 phase 2: 12-fragment evidence features and linear-part window choice for the single decoy-whitened network (written before any phase 2 training)

Context: `RESULT_V4P1.md` (no phase 1 arm passes; D11 stays the base) and `PLAN_V4.notes` P2.

- **Window choice:** the linear part of D11 alone picks the right window (against the DIA-NN RT, reference only) for 2255 / 2286 precursors. The full D11 score picks it for 1971 / 2017.
- **12-fragment arrays:** `build12.py` reads the DIALibGen 0.13.1 export of the existing windows. The 6 original fragments in it reproduce `xic.npy` (raw 01: median relative difference 0, p99 0.037).
- **Feature check (raw 01, before this plan, reference only):** `f12_libcos` picks the right window over the wrong winner in 0.867 of stage-C precursors. Its AUC for stage D against the top 1 % of decoy winners is 0.901.

## Features

`features12.py`, label-free, 6 columns:
- `f12_support`;
- `f12_coelution`;
- `f12_libcos`;
- `f12_extra_support`;
- `f12_extra_coelution`;
- `f12_extra_share`.

Decoys get their target's predicted intensities, with m/z computed on the decoy sequence. This is the rule of the original library.

## Arms

All arms use the D11 settings (PLAN_V3_P2b). The trainer is `train_v4.py` with α = 0, the same seeds, population, folds and rounds.

| Arm | Change from D11 |
|---|---|
| D11F | z = 11 features + 6 f12 features |
| D11W | the D11 model unchanged (weights reproduced, `diag_readout.py`); each precursor's window = argmax of the linear part z·w_z; its score = the full D11 score of that window |
| D11FW | D11F, with the D11W window rule |

**W rule.** Each window's score is computed from the one trained model. Within a precursor, the window with the highest linear part is kept, and its full score is the precursor score. The other windows get −∞. The rule is the same for targets and decoys and uses no labels. There is no second model and no extra training. The rounds' positive and negative selection inside training is unchanged; the W rule applies only to the final scores.

## Pass rule

Each arm must pass on both raws. The rules are the same as PLAN_V4.

1. PLAN_V2 rules 1–4 against BASE:
   - the ΔN@2 % interval is above 0;
   - the ΔN@1 % interval is not entirely below 0;
   - FDP_c at q ≤ 0.01 is at most 0.015;
   - the entrapment AUC is at most BASE's.
2. N@2 % is above D11 on both raws (paired bootstrap interval against D11 reported).
3. Null-winner AUC ≤ D11's + 0.005.

**Reported only:**
- the funnel against DIA-NN (stages C and D);
- N at 0.5 / 1 / 1.5 %;
- the window-choice rate.

An arm that passes becomes the base and the candidate for the C++ `OpenDIAScoringModel`.

## Deferred

Per-run fine-tuned RT / MS2 deviation features are not in this plan. DIALibGen RT matches the library RT (correlation 1.0), so it adds nothing, and a per-run fine-tune needs its own plan.
