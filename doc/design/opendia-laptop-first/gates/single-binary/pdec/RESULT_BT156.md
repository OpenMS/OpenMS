# RESULT_BT156: hard-decoy training holds at r ≈ 1.56 on the current default (21 Nov 2026, diagnostic)

The gate is in `PRECHECK_BT156.md` (`c52508fd28`), and the script is `precheck_ceil0.py` (sha256 `19955add…`) with r = 1.1225. The runs are in `$C/bt156-20261121`, current default, `-features_out`, and all exited 0. The full output is in `EVAL_BT156.txt`.

| Run | A tool N@FDP_c 1 % | Bt N@FDP_c 1 % | Gain | Lower bound A → Bt | HeLa at own q ≤ 0.01, A → Bt |
|---|---:|---:|---:|---|---|
| D s0 | 16,678 | 17,434 | +4.5 % | 0.89 → 0.94 % | 17,740 → 18,832 (+6.2 %) |
| D s1 | 17,038 | 17,524 | +2.9 % | 0.88 → 0.89 % | 17,929 → 18,713 (+4.4 %) |
| H s0 | 32,310 | 34,378 | +6.4 % | 0.66 → 0.75 % | 32,505 → 35,023 (+7.7 %) |
| H s1 | 32,415 | 34,244 | +5.6 % | 0.68 → 0.79 % | 32,622 → 35,116 (+7.6 %) |

| Criterion | D | H | |
|---|---|---|---|
| Mean gain at N@FDP_c 1 % ≥ +3 % | +3.7 % | +6.0 % | yes |
| Every run gains | 2 / 2 | 2 / 2 | yes |
| Mean lower-bound rise ≤ +0.2 pp | +0.03 pp | +0.10 pp | yes |

**Verdict: PASS.** Bt goes to a C++ opt-in port, then a two-seed check, then a frozen r7 plan.

- Plain label-free B is consistently weaker (D +3.5 % / +1.4 %, H +4.3 % / +3.9 % at N@FDP_c 1 %), so the hard-decoy negatives carry the gain.
- The DIA-NN-labelled ceiling C is below Bt on all four runs. This matches RESULT_CEIL0: the gain comes from sharpening the decoy tail, not from new information.
- **Why this differs from PRECHECK_NOCAND section 3** (+1.2 %, with the lower bound rising every round). That test ran on the older default (36 features, no matched filter, no RT refinement) and on one seed of D only. On the current 45 features the same design is stable across seeds and runs.

## The design that the C++ port must reproduce

`precheck_ceil0.py`, arm Bt:

1. Round 0 scores each window with the tool score.
2. Rounds 1 to 3, each cross-fitted over the library's 3 folds:
   - Run best-window-per-precursor paired TDC.
   - Positives are the target winners at q ≤ 0.01.
   - Negatives are all windows of the hardest 2 % of decoy precursors, ranked by their best-window score under the current round's score.
   - Fit a GBM (300 rounds, 31 leaves, balanced weights). The next round uses its out-of-fold scores.
3. Final q comes from paired TDC on the round-3 score.

The C++ `CrossFitScorer` is a small network with a decoy-whitened readout, not a GBM. Two things therefore have to be checked before the r7 plan: that the gain survives the change of learner, and that hard negatives fit the whitened-readout objective.
