# PLAN_GBT: hard-negative GBM rescoring on the r ≈ 7 entrapment (frozen 22 Nov 2026)

**Question.** Does `-gbm_rescore` (with `-gbm_hard_fraction 0.02`) raise IDs at equal true FDP on the large library without materially loosening the tool's q values? On top of the sibling default, it adds three cross-fitted histogram-GBM relabel rounds. Each round trains q ≤ 0.01 target winners against all windows of the hardest 2 % of decoy precursors.

## Evidence before this freeze

- **RESULT_GAP2.** On the sibling default at r7, seed 0, the Python Bt arm gains D +12.6 % and H +7.4 % at N @ FDP_c 1 %.
- **RESULT_HNEG.** The C++ port on r ≈ 1.56, seeds 0 and 1, gains:
  - **N @ FDP_c 1 %:** D +4.7 %, H +5.9 %, and every seed gains;
  - **lower bound:** up +0.10 pp on both runs;
  - **pair fairness:** D 0.542, H 0.518;
  - **wall:** D 1.49×, H 0.90×.

  With the flag off, output is byte-identical to the sibling default on all 4 runs.
- **Code.** OpenMS `d8ca86edeb`, 27/27 tests.
- **Disclosed:** these bars are the standard PLAN_SIB bars. They were set after RESULT_HNEG was seen, and they are not loosened for it. The D wall at r ≈ 1.56 (1.49×) and the D pair fairness (0.542) are both close to their bars.

## Setup (fixed)

- **Container:** `/nfs/research/juan/singj/opendia/opendia-pd10-gbmopt.sif`, sha256 `9b5df05f0972648ddd0c5ac7d9da2a9bf911c6b0ca810582935dd2fd5d2bada8`. Source is OpenMS `d8ca86edeb`.
- **Library:** `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`, as in PLAN_SIB.
- **Arms:**
  - **CUR:** the current default (siblings on). Its `precursors.tsv` must be byte-identical to `entlib-sib-20261121/runs/<R>_sib_s<seed>`, the adopted default. Otherwise the plan stops and nothing is evaluated.
  - **GBT:** CUR plus `-gbm_rescore -gbm_hard_fraction 0.02`. Both arms use the same container.
- **Runs:** D seeds 0-4, 8 threads, 128 GB. H (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py`, sha256 `bb03387dceccf01d4f446904cb8bf2f89deb9542fdcc21cdde503ad664b2b8a1`.

## Criteria (on D, pooled or median over the 5 seeds)

1. **Sensitivity:** GBT median N at FDP_c 1 % ≥ 1.03 × CUR median.
2. **Calibration not materially worse:**
   - GBT pooled lower bound at own q ≤ 0.01 ≤ CUR + 0.0020.
   - GBT pooled combined bound ≤ CUR + 0.0030.
3. **Pair fairness:** GBT P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.
4. **Cost:** GBT median wall ≤ 1.5 × CUR median wall, and peak RSS ≤ 64 GB.

No retuning, seed removal or post-hoc filtering.

## Outcome

- **All pass:** report H. Making `-gbm_rescore` the default is a separate commit with a byte check against the tested arm.
- **Criterion 1 fails:** the flag stays opt-in. Nothing is retuned.
- **Criterion 2, 3 or 4 fails with criterion 1 passing:** report the gain with its cost. The flag stays opt-in.
