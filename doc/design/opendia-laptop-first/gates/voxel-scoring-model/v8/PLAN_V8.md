# Voxel v8: the D11FWT recipe on the linear-end-RTMap population (written before any rtlin scores were computed)

## Context

- C++ commit `6015a6e53f` makes `OpenDIARTMap` extrapolate linearly beyond the end knots and is always on.
- The search-level test (`RTX_CHECK.md`) showed fewer A and B misses with linear ends.
- `rtlin-20261111/chain_data.sh` rebuilds everything with the new binaries (`binary.sha256`): calibration, search, seeds, BASE (OpenDIACNN crossfit `-stage2` search+boundary+byion), population, 12-fragment library and arrays.
- The question this plan answers: does the D11FWT recipe, unchanged, still pass against the rebuilt BASE on the new population, and how do N and the miss funnel compare with the old population?

## Arm (no new features or parameters; same scripts, paths repointed to `rtlin-20261111`)

D11FWT:
1. `train_v4.py <raw> D11F`, with `TRAIN_THREADS=10`;
2. `make_w.py <raw> D11F lin_D11F D11FW`;
3. `calib_tail.py <raw> D11FW D11FWT`.

Seeds, folds and hyperparameters are as in PLAN_V4 phase 2, PLAN_V6 and the scripts.

## Comparison base

**BASE** = the rebuilt CNN + LDA `stage2_score` from `rtlin-20261111/raw0?/byion/features.tsv`, restricted to the new population. Same metric and TDC (`evaluate_v2.py`).

## Pass rule

PLAN_V2 rules 1–4 against the rebuilt BASE, on both raws:
1. The bootstrap 90 % interval of N@2 % minus BASE lies above 0.
2. The bootstrap 90 % interval of N@1 % minus BASE does not lie entirely below 0.
3. FDP_c at q ≤ 0.01 is ≤ 0.015.
4. The entrapment-over-decoy winner AUC is ≤ the rebuilt BASE AUC.

**Reported only, descriptive:**
- N at 0.5 / 1 / 1.5 / 2 %, q ≤ 0.05;
- null-winner AUC;
- the A–E miss funnel against DIA-NN (`misses.py`);
- the old-population D11FWT numbers: N@1 % 1496 / 1426, N@2 % 1716 / 1736, FDP_c at q ≤ 0.01 0.0137 / 0.0131;
- old and new BASE.

**Caveat.** The populations differ, so old against new N is not a paired comparison.

**Outcome.**
- If it passes, D11FWT on the rtlin population becomes the base, and the C++ RTMap change stays.
- If it fails, the result is recorded and the next step is the user's decision. The RTMap change is not reverted without the user.
