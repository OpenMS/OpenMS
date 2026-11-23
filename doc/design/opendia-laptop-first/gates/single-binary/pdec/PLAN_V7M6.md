# PLAN_V7M6: candidate split `-valley_fraction 0.7 -max_peaks 6` on the r ≈ 7 entrapment (frozen 23 Nov 2026)

**Question.** Does V7M6 raise IDs at equal true FDP on the large library without materially loosening the tool's q values? V7M6 keeps two maxima apart up to a valley of 0.7 of the lower one (default 0.5) and keeps 6 candidate peaks per precursor (default 5), on top of the current default (siblings, GBM, faster interference resolver).

## Evidence before this freeze

- **RESULT_GAP2** (sibling default, seed 0): D r7 +13.3 %, H r7 +3.6 % at N @ FDP_c 1 %. It failed rule B on wall time only (1.52×).
- **RESULT_V7M6T** (GBM default, matched nodes): D r7 +13.3 %, H r7 +4.2 %, lower bound −0.05 pp. It failed rule B on D r7 wall (1.504×), and 62 % of the extra time was the interference step.
- **RESULT_V7M6T2** (faster resolver, OpenMS `955516d84c`, output byte-identical): the same gains, and rule B **passes** with D r7 wall 1.375×.
- **Disclosed:** all of this is seed 0. The bars below are the standard PLAN_SIB / PLAN_GBT bars, not set or loosened for V7M6.

## Setup (fixed)

- **Container:** `/nfs/research/juan/singj/opendia/opendia-pd10-fastint.sif`, sha256 `8de3a9684b1c0abc49ef0f385e7d9daca78e85f21184f51455c430633b9ee7d3`. Source is OpenMS `955516d84c`.
- **Library:** `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`, as in PLAN_SIB and PLAN_GBT.
- **Arms:**
  - **CUR:** the current default, no flags. Its `precursors.tsv` must be byte-identical to `entlib-gbt-20261122/runs/<R>_gbt_s<seed>`, the adopted GBM default. Otherwise the plan stops and nothing is evaluated.
  - **V7M6:** CUR plus `-valley_fraction 0.7 -max_peaks 6`. Both arms use the same container.
- **Runs:** D seeds 0-4, 8 threads, 160 GB, standard codon nodes (as PLAN_GBT). H (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py`, sha256 `bb03387dceccf01d4f446904cb8bf2f89deb9542fdcc21cdde503ad664b2b8a1`.

## Criteria (on D, pooled or median over the 5 seeds)

1. **Sensitivity:** V7M6 median N at FDP_c 1 % ≥ 1.03 × CUR median.
2. **Calibration not materially worse:**
   - V7M6 pooled lower bound at own q ≤ 0.01 ≤ CUR + 0.0020.
   - V7M6 pooled combined bound ≤ CUR + 0.0030.
3. **Pair fairness:** V7M6 P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.
4. **Cost:** V7M6 median wall ≤ 1.5 × CUR median wall, and peak RSS ≤ 64 GB.

No retuning, seed removal or post-hoc filtering.

## Outcome

- **All pass:** report H. Making V7M6 the default (changing the `-valley_fraction` and `-max_peaks` defaults) is a separate commit with a byte check against the tested arm.
- **Criterion 1 fails:** the defaults stay. Nothing is retuned.
- **Criterion 2, 3 or 4 fails with criterion 1 passing:** report the gain with its cost. The defaults stay.
