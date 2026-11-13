# RESULT_DECOYGEN: decoys with their own predicted spectrum and RT remove the entrapment excess (rule 1c passes in all arms)

- **Run:** PLAN_DECOYGEN, frozen 2026-11-13T20:26:54Z.
- **Amendment 1** (20:50Z, before any search output): the PAIRLIB run calibration is reused per arm.
- **Jobs:** codon 27900614–27900625 (`jobs_amend1.txt`).
- **Wall time:** about 55–60 min per arm, with the 4 arms in parallel. Search was about 35 s per raw. BASE crossfit took
  5:51 and 10 GB on raw01.
- **Decoy prediction:** DIALibGen on one codon A100, 35 s per prediction for 1.73M decoy precursors.
- **Outputs:** `results/arm*/`.

## Pass rule (D11FDWT; PLAN_XIC_IDS rules, unchanged)
| rule | A shuffle | B mutate | C pseudo-rev, no swap | D pseudo-rev + swap |
|---|---|---|---|---|
| N@2 % − BASE 90 % > 0 (raw01 / raw02) | pass / pass | pass / pass | pass / pass | pass / pass |
| N@1 % − BASE not all < 0 | pass / pass | pass / pass | pass (−379..+989) / pass | pass / pass |
| **1c FDP_c at q ≤ 0.01 ≤ 0.015** | **0.002 / 0.010 pass** | **0.004 / 0.003 pass** | **0.011 / 0.009 pass** | **0.008 / 0.008 pass** |
| AUC(entrap > decoy) ≤ BASE | 0.4980 ≤ 0.5031 / 0.5001 ≤ 0.5023 | **0.5014 > 0.4916 / 0.5021 > 0.4925** | 0.4997 ≤ 0.5005 / 0.4994 ≤ 0.5016 | **0.5013 > 0.4993** / 0.4975 ≤ 0.4996 |
| 2 N@FDP_c 2 % > 1716 / 1736 | 1917 / 1880 | 1860 / 1868 | 1949 / 1846 | 1864 / 1801 |
| 3 null-winner AUC ≤ 0.4983 / 0.4986 | 0.4970 / **0.4989** | 0.4940 / 0.4942 | 0.4974 / **0.4987** | 0.4972 / 0.4977 |

- **Verdict:** no arm passes every rule, so all four are formally FAIL. Rule 1c, the blocker of every earlier run,
  passes in all four arms on both raws.
- **Where the failures are:** the rule 3 misses (A, C) are 0.0001–0.0003. D's rule 1 AUC miss is 0.002.
- **Why those two rules now misfire:** their limits were fixed from runs with the old decoys. With calibrated decoys both
  AUCs sit at 0.50 ± noise. These rules were meant to catch a scorer that separates entrapments from decoys, and that no
  longer shows. Recorded, not retuned.

## Hypothesis (pooled rho over both raws, D11FDWT; `RHO_D11FDWT.txt`)
| arm | q ≤ 0.01 | q ≤ 0.05 | q ≤ 0.10 | q ≤ 0.20 | H |
|---|---|---|---|---|---|
| PAIRLIB (reference) | 2.76 [1.87, 4.12] | 1.59 [1.35, 1.87] | 1.42 [1.27, 1.59] | 1.33 [1.24, 1.43] | — |
| A shuffle | 0.86 [0.40, 1.74] | 0.90 [0.73, 1.10] | 1.15 [1.01, 1.30] | 1.00 [0.92, 1.09] | not supported (q ≤ 0.10 interval excludes 1), not rejected |
| B mutate | 0.43 [0.16, 1.03] | 0.38 [0.28, 0.53] | 0.59 [0.49, 0.70] | 0.61 [0.55, 0.68] | not supported (too conservative), not rejected |
| C pseudo-rev, no swap | 1.07 [0.60, 1.84] | 0.83 [0.67, 1.03] | 0.96 [0.84, 1.10] | 0.97 [0.89, 1.05] | **supported** |
| D pseudo-rev + swap | 0.93 [0.50, 1.68] | 0.96 [0.77, 1.18] | 1.08 [0.94, 1.23] | 1.10 [1.01, 1.19] | **supported** |

- **The prediction was wrong on two counts:**
  - D (the current sequence rule) is as well calibrated as C and A, not "between PAIRLIB and A";
  - B is not anti-conservative but strongly conservative.

## Label-free null check (`NULLRATES.txt`; decoy / target ratio in null windows, raw01 / raw02)
| arm | support ≥ 4 | libcos > 0.8 | eligible (0.92–1.08) |
|---|---|---|---|
| PAIRLIB | 0.854 / 0.830 | 0.803 / 0.771 | no |
| A | 1.023 / 1.005 | 0.949 / 0.935 | yes |
| B | **1.181 / 1.179** | 0.941 / 0.957 | **no** |
| C | 1.027 / 1.028 | 0.961 / 0.938 | yes |
| D | 1.016 / 0.997 | 0.971 / 0.956 | yes |

- **Labelled split** (`NULLRATES_labels.txt`, evaluation only): in A, C and D, the H, E, D(H) and D(E) classes now match
  noise at the same rate.
- **B:** its decoys match noise more often than targets (support ≥ 4: 6.5 % against 5.5 %), which fits its conservative
  rho. Why this happens is not tested here.

## Choice for C++ (plan rule, step 2)
- **Step 1 is empty:** no eligible, supported arm passes every pass rule.
- **Step 2:** of the supported and eligible arms (C, D), D has the pooled q ≤ 0.05 rho closest to 1 (0.96 against 0.83;
  not a tie). **D is chosen:** pseudo-reverse with the K/R swap, but the decoy gets its own predicted fragment
  intensities, fragment selection and RT.
- **Step 2's note on rule 1c** does not apply: rule 1c passes. The remaining failures are the two near-0.5 AUC rules.

## Reading
- **The excess came from how the decoy spectrum and RT were assigned, not from the sequence rule.**
  - PAIRLIB and D share the same decoy sequences (98.8 % identical).
  - They differ only in the decoy's intensities, fragment choice and RT: copied from the target in PAIRLIB, predicted
    from the decoy's own sequence in D.
  - That alone moves the null-match ratio from 0.83 to 1.0, pooled rho at q ≤ 0.05 from 1.59 to 0.96, and FDP_c at
    q ≤ 0.01 from 0.032 / 0.017 to 0.008 / 0.008.
- **Mechanism, not tested:** copied intensities in a reversed ladder put the target's high-intensity slots on the decoy's
  ion types and ordinals that are rarely intense. So the decoy's top-6 m/z set is less like a real peptide's, and it
  hits real background less often.
- **Not separated here:** the predicted intensities against the predicted RT. D changes both.
- **IDs:** N@FDP_c 1 % for D11FDWT is 1733 / 1618 for D, against 1549 / 1671 for PAIRLIB.
  - At the tool's own q ≤ 0.01, PAIRLIB reported 1978 / 1784 at FDP_c 0.032 / 0.017; D reports 1667 / 1445 at
    0.008 / 0.008.
  - So reported IDs are now honest, and the ranking gap to DIA-NN remains: DIA-NN 2.2 on PAIRLIB raw01 gives 2067 at
    N@FDP_c 1 %.
- **Finding (amendment 1):** OpenDIACalibrate failed or degraded with the new decoy sets, from 607 anchors down to
  0–515. The calibrator depends on its decoy set more than it should.

## For the C++ implementation
- The generator only needs to predict decoys from their own sequence.
  - OpenSwathDecoyGenerator and DIALibGen's `-generation:decoys` both keep the target's intensities. DIALibGen records
    this as `decoy_semantics: target-sequence-with-shifted-fragments`, which is the old behaviour.
  - The decoy needs its own MS2 and RT prediction, and its own top-N fragment selection, at the target's charge.
- The DIALibGen prediction path already does this when given the decoy sequence as a peptide, which is how this run
  built the libraries. A predicted-decoy mode is mostly plumbing.
- The label-free null-rate check (`null_rates.py`) is the deployable acceptance test for any decoy generator.
