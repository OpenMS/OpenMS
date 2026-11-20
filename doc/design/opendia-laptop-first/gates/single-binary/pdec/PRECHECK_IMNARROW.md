# PRECHECK_IMNARROW: a narrower calibrated 1/K0 window (20 Nov 2026, diagnostic)

This is a diagnostic. It led to PLAN_IMNARROW, which was frozen before any r ≈ 7 run. DIA-NN and entrapment labels
are used for evaluation only.

**Idea.** A faint peak's fragment traces collect every point inside the calibrated mobility window, ±0.052 1/K0
(4 σ of the anchor residual). Most of that window is background for a faint peak. Narrowing it should raise signal to
background for real peaks. A decoy gains nothing, since it has no consistent mobility, and it goes through the same code.

**Change.** A new advanced option, `-im_tolerance_scale` (0942ffbdf3), multiplies the calibrated half width. The
scaled width is used by the search, the 12-fragment grid, the interference resolver and the evidence features. The
6-fragment trace tensor keeps its fixed ±0.03. At 1 the output is byte-identical to the previous default on
r ≈ 1.56 HeLa D, and 25/25 ctests pass.

## Sweep on r ≈ 1.56 HeLa D

Seed 0 unless marked. `eval_pdec.py` evaluates the tool output as written. Scripts: `imnarrow.sh`, `imnarrow2.sh`,
`precheck_imnarrow.py`, `precheck_arms.py`.

| scale (half width) | N at FDP_c 1 % | q ≤ 0.01 lower / combined bound | pair win | wall / RSS |
|---|---|---|---|---|
| 1.0 (±0.052, default) | 12,357 | 0.79 / 1.49 % | 0.524 | 5:03 / 11.0 GB |
| 0.8 | 12,806 (+3.6 %) | 0.69 / 1.31 % | 0.507 | 5:08 / 11.1 GB |
| **0.6 (±0.031)** | **13,388 (+8.3 %)** | **0.57 / 1.08 %** | 0.501 | 6:13 / 11.2 GB |
| 0.5 | 13,666 (+10.6 %) | 0.60 / 1.13 % | 0.518 | 6:11 / 11.8 GB |
| 0.4 | 13,103 (+6.0 %) | 0.72 / 1.37 % | 0.510 | 7:38 / 12.4 GB |
| seed 1: 1.0 → 0.6 (check only, after the freeze) | 12,339 → 13,545 (+9.8 %) | 0.74 / 1.40 → 0.64 / 1.22 % | 0.506 → 0.517 | 5:34 → 7:02 |

- More IDs and better calibration together. Fewer entrapments are accepted (78 against 103 at q ≤ 0.01) while N rises. So the gain does not come from loosening.
- Below 0.5 the window starts to cut real signal, since the true error has p95 0.031, and the gain falls.
- The extra cost is mostly more candidate windows (2.39 M to 2.66 M at 0.6), because narrower traces split fewer peaks together.

**Faint peaks gain most.** Share of all DIA-NN IDs whose right window, scored by the tool, beats the decoy tail:

| Decoy FPR | Arm | Q1 | Q2 | Q3 | Q4 | all |
|---|---|---|---|---|---|---|
| 1e-4 | 1.0 | 0.143 | 0.262 | 0.472 | 0.801 | 0.420 |
| 1e-4 | 0.6 | 0.166 | 0.298 | 0.504 | 0.809 | 0.444 |
| 1e-2 | 1.0 | 0.312 | 0.451 | 0.637 | 0.887 | 0.572 |
| 1e-2 | 0.6 | 0.352 | 0.505 | 0.681 | 0.893 | 0.608 |

- A right candidate is present for 84.1 % of DIA-NN IDs at 1.0 and 89.1 % at 0.6. This is the same masking that PRECHECK_IMPICK found, removed here without a 2-D picker.
- IDs accepted with the right window: 47.3 % → 49.7 % (1,013 gained, 402 lost).

## Decision

This is the largest single lever since PLAN_EXTRA. **PLAN_IMNARROW** was frozen (07a4245472) with scale 0.6 for the r ≈ 7 D 5-seed test. That plan states how 0.6 was chosen and the weaknesses of that choice.
