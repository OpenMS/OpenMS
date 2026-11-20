# PLAN_MT: fragment-level matched-filter features on the r ≈ 7 entrapment (frozen 20 Nov 2026)

**Question.** Does `-matched_filter`, 9 per-fragment robust features computed from the existing trace tensor, raise IDs at equal true FDP on the large library without materially loosening the tool's q values?

## Evidence before this freeze

- **PRECHECK_MF** (docs 7b7a86f842). This was a raw-point precheck on the faint2 sample, with its gate written before any result. The tensor variant `mt` raised the label-free proxy TPR at decoy FPR 3e-3 on missed DIA-NN IDs by +0.103 / +0.091 over the 36 default features, against +0.025 / +0.013 for core5.
- **C++** is OpenMS `34c36861f2`. It is opt-in, and 25/25 OpenDIA tests pass. With the flag off, the output is byte-identical to the current default on r ≈ 1.56 D, seeds 0 and 1: `precursors.tsv`, the run table and `calibration.json` all match the `faintcpp-20261120` default runs.
- **Workflow check, r ≈ 1.56 D**, one binary, `eval_pdec.py`:

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % |
|---|---|---:|---:|---:|---:|---:|
| CUR | 0 | 15,899 | 130 | 0.82 % | 1.55 % | 14,948 |
| CUR | 1 | 15,798 | 131 | 0.83 % | 1.57 % | 14,867 |
| MT | 0 | 17,473 | 149 | 0.85 % | 1.61 % | 16,385 |
| MT | 1 | 17,221 | 141 | 0.82 % | 1.55 % | 16,101 |

  - The median N at FDP_c 1 % goes from 14,908 to 16,243 (+9.0 %), and every MT seed is above every CUR seed.
  - Pooled bounds go from 0.82 / 1.56 % to 0.84 / 1.58 %, and pair win from 0.523 to 0.514.
  - Wall time is 2:39 against 2:35 for the same-node seed 0, and peak RSS is about 9.3-9.5 GB for both arms.
- **Disclosed:** PLAN_FAINT shrank from +6.3 % at r ≈ 1.56 to +4.0 % at r ≈ 7. The bar below is the same +5 % as PLAN_FAINT and was not tuned to this result.

## Setup (fixed)

- **Container:** `/nfs/research/juan/singj/opendia/opendia-pd10-34c36861f2.sif`, sha256 `1cad85f5605b3c4a62ee7d860594aaaddcdcd1dc9099bd930ad5446f1ad16d06`. Source is OpenMS `34c36861f2`.
- **Library:** the PLAN_ODP predicted r ≈ 7 library, `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`, as in PLAN_FAINT.
- **Arms:** CUR is the current default. MT is CUR plus `-matched_filter`. Both arms use the same container.
- **Runs:** D seeds 0-4, 8 threads, 128 GB. H (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py`, sha256 `bb03387dceccf01d4f446904cb8bf2f89deb9542fdcc21cdde503ad664b2b8a1`. Entrapment labels are read only there.

## Criteria (on D, pooled or median over the 5 seeds)

1. **Sensitivity:** MT median N at FDP_c 1 % ≥ 1.05 × CUR median.
2. **Calibration not materially worse:**
   - MT pooled lower bound at own q ≤ 0.01 ≤ CUR + 0.0020.
   - MT pooled combined bound ≤ CUR + 0.0030.
3. **Pair fairness:** MT P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.
4. **Cost:** MT median wall ≤ 1.3 × CUR median wall, and peak RSS ≤ 64 GB.

No retuning, feature removal, seed removal or post-hoc filtering.

## Outcome

- **All pass:** report H. Making `-matched_filter` the default is a separate commit with ctests and a byte check against the tested arm.
- **Criterion 1 passes and criterion 2 fails:** report the gain with its calibration cost. The flag stays opt-in, and nothing is retuned.
- **Any other fail:** report as is. The flag stays opt-in.
- **Follow-up either way:** the raw-point `mf` variant had about twice the ceiling gain of `mt` in PRECHECK_MF. It needs a finer, ppm-weighted grid in the extraction pass and gets its own plan.
