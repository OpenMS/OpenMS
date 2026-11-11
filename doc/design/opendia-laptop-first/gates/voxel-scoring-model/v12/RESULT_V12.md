# V12 result: interference-corrected profile view

## Decision

**FAIL. Keep D11FWT as the base. Do not promote D11FIWT.** The two-view voxel arm passes the null-winner guard, but misses the frozen FDR requirement on raw01 and returns fewer identifications at the 2% count point than D11FWT on both raws.

## Frozen protocol and execution

- Arm: `D11FIWT` (D11FWT inputs plus one dual-view voxel scorer).
- Comparison base: `D11FWT`; raws 01 and 02 only.
- Frozen protocol: `PLAN_V12.md`, SHA-256 `58197a68b92358c3a5097f968cdbf13771f2e81394a5f02bd3e7cefa28a7a196` (written at 2026-11-11T22:12:29Z).
- Frozen source checks: `PLAN_V12.code.sha256`; the chain verified the plan before running.
- Chain status: `done`, exit code `0`.
- Training was CPU-only. Each raw's three-fold `D11FI` fit took about 35 minutes (raw01 34:55.63, raw02 35:19.91), followed by W, tail calibration, and evaluation. The first raw01 attempt was deliberately stopped by PID while checking whether the frozen run could move to codon; it was not resumed. The chain was then restarted from the frozen plan and completed successfully. No model changes were made after the plan freeze.

## Frozen pass gates

1. **PLAN_V2 rules 1–4 against BASE, including FDP_c at q ≤ 0.01 ≤ 0.015: FAIL.** Raw01 D11FIWT FDP_c is 0.027. Raw02 is 0.010.
2. **N@2% greater than D11FWT on both raws: FAIL.** Raw01 is 1,585 vs 1,716; raw02 is 1,623 vs 1,736.
3. **Null-window winner AUC no more than D11FWT + 0.005: PASS.** Raw01 is 0.4939 vs 0.4933 (+0.0006); raw02 is 0.4930 vs 0.4936 (-0.0006).

The frozen combined rule therefore fails.

## Evaluation

| Raw | Arm | N@0.5% | N@1% | N@1.5% | N@2% | q ≤ 0.01 N | FDP_c at q ≤ 0.01 | Entrapment-vs-decoy AUC |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 01 | BASE | 900 | 1,024 | 1,267 | 1,314 | 1,081 | 0.013 | 0.5313 |
| 01 | D11FWT | 971 | 1,496 | 1,616 | 1,716 | 1,614 | 0.014 | 0.4907 |
| 01 | D11FIWT | 1,076 | 1,433 | 1,544 | 1,585 | 1,615 | 0.027 | 0.4878 |
| 02 | BASE | 1,215 | 1,288 | 1,333 | 1,376 | 1,241 | 0.008 | 0.5348 |
| 02 | D11FWT | 1,127 | 1,426 | 1,625 | 1,736 | 1,536 | 0.013 | 0.4958 |
| 02 | D11FIWT | 955 | 1,466 | 1,576 | 1,623 | 1,445 | 0.010 | 0.4888 |

The 90% paired-bootstrap intervals for D11FIWT versus BASE at the 1% point are raw01 [-2, +604] and raw02 [-304, +341]. At 2%, they are raw01 [+196, +446] and raw02 [+135, +315]. Versus D11FWT, the 90% intervals at 1% are raw01 [-396, +449] and raw02 [-355, +153]; at 2% they are raw01 [-189, -9] and raw02 [-199, +47]. These intervals do not change the frozen decision because the arm fails the explicit FDP and N@2% gates.

## Interpretation

The model could consume both raw and corrected chromatograms, but the learned combination did not improve the target outcome. On raw01 the arm reached nearly the same q≤0.01 count as D11FWT while exceeding the allowed FDP_c; on both raws its N@2% was lower. Null AUC stayed within the allowed margin, so this is not a null-separation regression under the frozen check, but that alone is insufficient to accept the arm.

This is evidence against this particular reference-trace cap plus dual-view model experiment, not proof that chromatographic interference is irrelevant or that a voxel model cannot learn it. The pre-plan feature differences were modest, and the correction can remove legitimate signal. Stop this correction approach as planned and retain D11FWT. Do not alter the comparison base, pass rule, or promote C++ behavior based on this failed arm.

## Artifacts

- `v4/EVALUATION_V12.txt`
- `v4/NULL_V12.txt`
- `v4/VS_D11FWT_V12.txt`
- `v4/funnel_D11FIWT.txt`
- `v4/train_v12_raw01.log`, `v4/train_v12_raw02.log`
- `v4/make_w_v12.log`, `v4/calib_tail_v12.log`
