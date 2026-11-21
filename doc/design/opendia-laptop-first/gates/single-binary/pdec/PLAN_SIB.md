# PLAN_SIB: charge-state sibling rescoring on the r ≈ 7 entrapment (frozen 21 Nov 2026)

**Question.** Does `-sibling_features` raise IDs at equal true FDP on the large library without materially loosening the tool's q values? The flag scores once, adds `s_n`, `s_co`, `s_best` and `s_best_d` per window from those calibrated out-of-fold scores, then scores again.

## Evidence before this freeze

- **PRECHECK_JOINT / JOINT2.** Siblings carry new information. The GBM ceiling rises by +21.5 % on D and +9.8 % on H at r0.
- **RESULT_SIBCPP.** On r ≈ 1.56, seeds 0 and 1, the C++ flag gains D +38.9 % and H +13.5 % at N @ FDP_c 1 %. Every seed gains, and the lower bound changes by −0.13 pp on D and −0.01 pp on H.
- **Default unchanged.** With the flag off, output is byte-identical on all 4 runs.
- **Fairness and cost.** Pair fairness goes from 0.523 to 0.499 on D and from 0.511 to 0.507 on H. Wall time and RSS are within node noise.
- **Code.** OpenMS `a3eea6d9ee`, 26/26 tests.
- **Disclosed:** the +3 % bar below is the standard bar from earlier plans and was set before RESULT_SIBCPP. At r ≈ 7 the first-pass scores, and so `s_best`, are noisier, so the gain may shrink (the matched filter went from +9.0 % to +6.5 %).

## Setup (fixed)

- **Container:** `/nfs/research/juan/singj/opendia/opendia-pd10-sibf.sif`, sha256 `c67d7d9fa8cf0f127f8326357ee02ec8d34b14d6922e881250331ca8a23e28f6`. Source is OpenMS `a3eea6d9ee`.
- **Library:** `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`, as in PLAN_RTR.
- **Arms:**
  - **CUR:** the current default. Its `precursors.tsv` must be byte-identical to `entlib-rtr-20261121/runs/<R>_rtr_s<seed>`, the adopted default. Otherwise the plan stops and nothing is evaluated.
  - **SIB:** CUR plus `-sibling_features`. Both arms use the same container.
- **Runs:** D seeds 0-4, 8 threads, 128 GB. H (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py`, sha256 `bb03387dceccf01d4f446904cb8bf2f89deb9542fdcc21cdde503ad664b2b8a1`.

## Criteria (on D, pooled or median over the 5 seeds)

1. **Sensitivity:** SIB median N at FDP_c 1 % ≥ 1.03 × CUR median.
2. **Calibration not materially worse:**
   - SIB pooled lower bound at own q ≤ 0.01 ≤ CUR + 0.0020.
   - SIB pooled combined bound ≤ CUR + 0.0030.
3. **Pair fairness:** SIB P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.
4. **Cost:** SIB median wall ≤ 1.5 × CUR median wall, and peak RSS ≤ 64 GB.

No retuning, seed removal or post-hoc filtering.

## Outcome

- **All pass:** report H. Making `-sibling_features` the default is a separate commit with a byte check against the tested arm.
- **Criterion 1 fails:** the flag stays opt-in. Nothing is retuned.
- **Criterion 2, 3 or 4 fails with criterion 1 passing:** report the gain with its cost. The flag stays opt-in.
