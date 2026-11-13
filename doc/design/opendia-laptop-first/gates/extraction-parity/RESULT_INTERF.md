# RESULT_INTERF: fragment-deficit features (PLAN_INTERF), FAIL on rule 1c (raw01)

- **Run:** `interf-20261113/`, with `PLAN_INTERF.md` written at 2026-11-13T09:34:49Z.
  - The plan and script hashes were checked at chain start, and the chain was run once.
  - Artefacts are in `interf/`: plan, chain, scripts and score outputs.
- **Arm D11FDWT:** D11F plus 3 label-free fragment-deficit columns (`d_max`, `d6_max`, `d_lib` from
  `deficit_features.py`), then W and T as for D11FWT.
- **Control:** the control `D11FWT_rerun` reproduces D11FWT exactly on both raws, because training is deterministic.
  Every difference below is therefore caused by the 3 columns.

## Rules (PLAN_XIC_IDS rules, unchanged)
| rule | raw01 | raw02 |
|---|---|---|
| 1 bootstrap 90 % of N@2 % − BASE > 0 | +576..+934 pass | +854..+1259 pass |
| 1 N@1 % − BASE not entirely < 0 | +544..+946 pass | +671..+1302 pass |
| **1c FDP_c at q ≤ 0.01 ≤ 0.015** | **0.023 FAIL** (N 2007, E 23) | 0.008 pass (N 1668, E 7) |
| 1 AUC(entrap > decoy winner) ≤ BASE | 0.4928 ≤ 0.5190 pass | 0.4979 ≤ 0.5152 pass |
| 2 N@FDP_c 2 % > 1716 / 1736 | 1954 pass | 1969 pass |
| 3 null-winner AUC ≤ 0.4983 / 0.4986 | 0.4971 pass | 0.4972 pass |

**Verdict:** FAIL. Rule 1c fails on raw01. By plan, there is no retuning.

## Recorded, not ruled on
**N@FDP_c 0.5/1/1.5/2 %**
- This measure is calibrated by entrapments, so it does not depend on how good the decoy model is.

| arm | raw01 | raw02 |
|---|---|---|
| D11FWT | 1414 / 1567 / 1716 / 1862 | 1302 / 1558 / 1651 / 1751 |
| D11FDWT | 1640 / 1800 / 1896 / 1954 | 1429 / 1764 / 1883 / 1969 |
| gain | +226 / +233 / +180 / +92 | +127 / +206 / +232 / +218 |

**E, D and rho = E(1+1/r)/D among winners** (`DIAG_D11FDWT.txt`)

| q cut | raw01 E / D / rho | raw02 E / D / rho |
|---|---|---|
| 0.01 | 23 / 19 / 2.44 [1.40, 4.28] | 7 / 15 / 0.94 [0.38, 2.14] |
| 0.05 | 109 / 122 / 1.80 | 93 / 116 / 1.62 |
| 0.10 | 217 / 280 / 1.56 | 194 / 264 / 1.48 |
| 0.20 | 488 / 691 / 1.42 | 487 / 680 / 1.44 |

**Entrapment winners at q ≤ 0.01 with a co-eluting same-composition reference peptide** (`PARENT_D11FDWT.txt`)
- raw01: 7/23 (D11FWT: 5/14).
- raw02: 3/7 (D11FWT: 4/14).

## Reading
1. **The deficit columns are real ranking evidence.**
   - At every fixed entrapment FDP, IDs rise on both raws, by +127 to +233. Raw02 is the replication, and it has
     the larger gain at 1.5 % and 2 %.
   - The null-window AUC and AUC(entrap > decoy) do not get worse.
2. **They do not fix the decoy calibration.**
   - rho at q ≤ 0.05 to 0.20 is 1.4 to 1.8, the same band as D11FWT and the old base in `RESULT_QCAL.md`.
   - The q ≤ 0.01 values (2.44 and 0.94) are counts of 7 to 23, and their intervals overlap that band.
   - Taking both raws together at q ≤ 0.01, FDP_c is 30·(1+1/r)/3675 = 0.016, which is the D11FWT value.
   - So rule 1c moved by count noise: up on raw01 and down on raw02. The excess of entrapments over decoys did
     not change.
3. **The parent-explained share did not drop.**
   - Among winners, the columns do not separate entrapments sitting on a parent peptide any better than D11FWT does.
   - The pre-plan raw01 AUC of 0.74 (score-matched, among all windows) did not carry over to the winners near the
     cut.
4. **Consequence for option A.**
   - A label-free scoring feature raises the score of true targets. It does not, by itself, make decoys sample the
     null that FDRBench shuffles occupy, so rho stays about 1.5.
   - Meeting rule 1c on both raws needs one of two things. One is a change to the decoy-to-entrapment correspondence,
     which is a population or decoy change and needs approval. The other is a target that is not dominated by
     count noise at q ≤ 0.01, which is a rule change.
   - Any new scoring plan with the same decoys is expected to land on FDP_c of about 1.5 × 0.01 again.
