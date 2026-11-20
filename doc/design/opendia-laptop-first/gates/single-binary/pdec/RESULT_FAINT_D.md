# RESULT_FAINT_D: PLAN_FAINT D verdict (20 Nov 2026)

## Verdict

**FAIL, criterion 1 only.** The faint trace evidence block increased the median N at FDPc 1 % by 4.0 %, below the frozen 5 % requirement. The calibration and cost safeguards passed.

No H run was submitted because the frozen D verdict does not pass sensitivity. The `-faint_features` block remains opt-in.

## Evaluation

The frozen evaluator is `pdec/eval_pdec.py`, sha256 `bb03387dceccf01d4f446904cb8bf2f89deb9542fdcc21cdde503ad664b2b8a1`. The r ≈ 7 D runs used the frozen container and seeds 0 through 4 for both arms.

| Arm | N q ≤ 0.01 | E q ≤ 0.01 | Lower | Combined | N at FDPc 1 % |
|---|---:|---:|---:|---:|---:|
| CUR seed 0 | 14,562 | 157 | 1.08 % | 1.23 % | 13,924 |
| CUR seed 1 | 14,670 | 162 | 1.10 % | 1.26 % | 14,267 |
| CUR seed 2 | 14,612 | 171 | 1.17 % | 1.34 % | 14,009 |
| CUR seed 3 | 14,578 | 159 | 1.09 % | 1.25 % | 13,985 |
| CUR seed 4 | 14,390 | 155 | 1.08 % | 1.23 % | 13,889 |
| FAINT seed 0 | 15,027 | 167 | 1.11 % | 1.27 % | 14,549 |
| FAINT seed 1 | 15,069 | 163 | 1.08 % | 1.24 % | 14,501 |
| FAINT seed 2 | 15,075 | 170 | 1.13 % | 1.29 % | 14,580 |
| FAINT seed 3 | 15,261 | 186 | 1.22 % | 1.39 % | 14,645 |
| FAINT seed 4 | 15,123 | 169 | 1.12 % | 1.28 % | 14,510 |

Pooled and median summaries:

| Arm | Pooled lower | Pooled combined | Median N at FDPc 1 % | Pair win |
|---|---:|---:|---:|---:|
| CUR | 1.10 % | 1.26 % | 13,985 | 0.499 |
| FAINT | 1.13 % | 1.29 % | 14,549 | 0.512 |

The sensitivity ratio is `14,549 / 13,985 = 1.0403`, versus the frozen minimum of 1.05. The lower-bound increase is 0.03 percentage points and the combined-bound increase is 0.03 percentage points, both within the frozen limits of 0.20 and 0.30 percentage points. Pair fairness passes. Same-node Slurm elapsed times had medians of 578 seconds for CUR and 560 seconds for FAINT. Peak RSS was approximately 33.6 GiB for CUR and 33.8 GiB for FAINT.

## Reading

The r ≈ 1.56 precheck predicted a 6.3 % median gain, but the five-seed r ≈ 7 result settled at 4.0 %. The direction is consistent and the cost is small, but it is not enough for the precommitted promotion rule. The result is retained as a useful opt-in diagnostic, not presented as a default improvement.
