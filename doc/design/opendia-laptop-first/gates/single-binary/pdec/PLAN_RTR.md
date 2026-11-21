# PLAN_RTR: per-run ridge RT refinement on the r ≈ 7 entrapment (frozen 21 Nov 2026)

**Question.** Does `-rt_refine` raise IDs at equal true FDP on the large library without materially loosening the tool's q values? It fits a Sage-style per-run ridge regression of the RT residual on sequence composition from the run's own q ≤ 0.01 IDs, then searches, extracts and scores again with corrected RTs and a narrower RT window.

## Evidence before this freeze

- **PRECHECK_RTIM** (docs 8929ffefbe). On missed DIA-NN IDs, the cross-fitted ridge correction lowers the p95 RT residual by 27 % on D and 38 % on H, below the best monotone map. A 1/K0 correction failed (−5 to −8 %), and the IM line is already at its isotonic floor (D 0.0310 against 0.0303, H 0.0302 against 0.0299), so mobility is not changed.
- **C++** is OpenMS `f115760e89`, which adds `-rt_refine`, plus `43da5895b1`, which adds `-rt_refine_anchors`. Both are opt-in, and 26/26 tests pass. With the flags off, the output is byte-identical to the current default (01c5891ce0) on r ≈ 1.56 D, seeds 0 and 1, in `precursors.tsv`, the run table, `calibration.json` and `features.tsv`.
- **Workflow check, r ≈ 1.56 D**, using `eval_pdec.py`:

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % |
|---|---|---:|---:|---:|---:|---:|
| default | 0 | 17,473 | 149 | 0.85 % | 1.61 % | 16,385 |
| default | 1 | 17,221 | 141 | 0.82 % | 1.55 % | 16,101 |
| `-rt_refine` | 0 | 17,839 | 154 | 0.86 % | 1.63 % | 16,678 |
| `-rt_refine` | 1 | 18,083 | 158 | 0.87 % | 1.65 % | 17,038 |
| `-rt_refine_anchors` | 0 | 17,516 | 159 | 0.91 % | 1.72 % | 16,491 |
| `-rt_refine_anchors` | 1 | 17,657 | 155 | 0.88 % | 1.66 % | 16,581 |

  - **`-rt_refine` median N at FDP_c 1 %:** 16,243 → 16,858, which is +3.8 %. Both seeds are above both default seeds.
  - **`-rt_refine` bounds and fairness:** pooled bounds go from 0.84 / 1.58 % to 0.87 / 1.64 %, and pair win from 0.514 to 0.523.
  - **`-rt_refine` refit:** it trained on about 17,300 IDs. The p95 residual dropped from 20.8-21.1 s to 15.6-15.7 s, which narrowed the RT window from ±35.2 s to ±26.2-26.4 s and cut candidates by 12 %.
  - **`-rt_refine` cost:** wall time is 1.7-1.9×, from the second pass. Peak RSS is about the same (9.5-10.1 GB).
  - **`-rt_refine_anchors`** trained on only 944 anchors and narrowed the window to ±29.9 s. It gained +1.9 % with a higher lower bound (0.91 / 0.88 %), so it is not taken.
- **Disclosed:** at r ≈ 1.56 the gain is +3.8 %, smaller than the matched filter's +9.0 %, which shrank to +6.5 % at r ≈ 7. The r ≈ 7 gain may be smaller still. The +3 % bar below was set after seeing this check.

## Setup (fixed)

- **Container:** `/nfs/research/juan/singj/opendia/opendia-pd10-rtrefine2.sif`, sha256 `5e532efea5ed227b7adbe85a16723d2b99ec987d701685fa316d355e7f18e38f`. Source is OpenMS `43da5895b1`.
- **Library:** the PLAN_ODP predicted r ≈ 7 library (`entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`), as PLAN_MT.
- **Arms:** CUR is the current default. RTR is CUR plus `-rt_refine`. Both arms use the same container.
- **Runs:** D seeds 0-4, 8 threads, 128 GB. H (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py`, sha256 `bb03387dceccf01d4f446904cb8bf2f89deb9542fdcc21cdde503ad664b2b8a1`.

## Criteria (on D, pooled or median over the 5 seeds)

1. **Sensitivity:** RTR median N at FDP_c 1 % ≥ 1.03 × CUR median.
2. **Calibration not materially worse:**
   - RTR pooled lower bound at own q ≤ 0.01 ≤ CUR + 0.0020.
   - RTR pooled combined bound ≤ CUR + 0.0030.
3. **Pair fairness:** RTR P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.
4. **Cost:** RTR median wall ≤ 2.2 × CUR median wall (the second pass is a full search and extraction), and peak RSS ≤ 64 GB.

No retuning, seed removal or post-hoc filtering.

## Outcome

- **All pass:** report H. Making `-rt_refine` the default is a separate commit with a byte check against the tested arm.
- **Criterion 1 fails:** the flag stays opt-in. The RT accuracy gain is reported as is, and nothing is retuned.
- **Criterion 2 or 3 fails with criterion 1 passing:** report the gain with its calibration cost. The flag stays opt-in.
