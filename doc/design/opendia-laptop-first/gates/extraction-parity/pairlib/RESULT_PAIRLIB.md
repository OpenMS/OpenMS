# RESULT_PAIRLIB: a library that contains every entrapment's parent does not remove the excess (FAIL, H rejected)

- **Run:** PLAN_PAIRLIB, frozen 2026-11-13T14:20Z, on codon (data jobs 27595228 and 27595229, scoring job 27630700).
- **Amendments:** two result-neutral ones.
  - Amendment 1: bind the pilot import files.
  - Amendment 2: give `parent.py` an absolute root. The diagnosis then runs locally on the copied outputs.
- **Wall time:** about 50 min on codon for both raws' data (in parallel), and 5 min for scoring. The same steps take
  about 3 h on the laptop.

## Pass rule (D11FDWT, PLAN_XIC_IDS rules)
| rule | raw01 | raw02 |
|---|---|---|
| 1 N@2 % − BASE 90 % above 0 | +502..+829 pass | +598..+1093 pass |
| 1 N@1 % − BASE not all below 0 | +85..+784 pass | +418..+990 pass |
| **1c FDP_c at q ≤ 0.01 ≤ 0.015** | **0.032 FAIL** (N 1978, E 32) | **0.017 FAIL** (N 1784, E 15) |
| 1 AUC(entrap > decoy) ≤ BASE | 0.4979 ≤ 0.5146 pass | 0.4946 ≤ 0.5193 pass |
| 2 N@FDP_c 2 % > 1716 / 1736 | 1821 pass | 1870 pass |
| 3 null-winner AUC ≤ 0.4983 / 0.4986 | 0.4964 pass | 0.4982 pass |

**Verdict: FAIL on rule 1c, both raws.**

## Hypothesis (ruled): pooled rho at q ≤ 0.05 ≤ 1.2 if missing parents caused the excess
| arm (pooled raw01 + raw02) | q ≤ 0.01 | q ≤ 0.05 | q ≤ 0.10 | q ≤ 0.20 |
|---|---|---|---|---|
| D11FDWT rho [90 %] | 2.76 [1.87, 4.12] | **1.59 [1.35, 1.87]** | 1.42 [1.27, 1.59] | 1.33 [1.24, 1.43] |
| D11FWT rho | 1.94 [1.24, 3.02] | 1.77 [1.50, 2.09] | 1.53 [1.36, 1.71] | 1.39 [1.29, 1.49] |

- rho at q ≤ 0.05 is 1.59, above the plan's rejection bound of 1.3. **H is rejected.**
- Putting every parent into the library leaves the excess where it was: 1.4–1.9 in earlier arms, 1.3–1.8 here.
- Missing parents are therefore not the main cause.

## Recorded, not ruled on
- **N@FDP_c 0.5/1/1.5/2 %** (entrapment-calibrated):

  | arm | raw01 | raw02 |
  |---|---|---|
  | D11FDWT | 1114 / 1549 / 1705 / 1821 | 1276 / 1671 / 1700 / 1870 |
  | D11FWT | 1239 / 1440 / 1553 / 1711 | 1302 / 1622 / 1754 / 1824 |

- **Same-composition claimants of entrapment winners** (raw01 D11FDWT; `interf_pre/COMPONENTS_pairlib_raw01.txt`):
  - 22 % at q ≤ 0.01 and 23 % at q ≤ 0.05 share signal with a candidate peak of the same composition, which is now
    their in-library parent.
  - Excluding those, raw01 rho at q ≤ 0.05 would still be about 1.5. Parent interference explains about a fifth of
    the excess, not most of it.
- **Shared signal by class** (≥ 2 shared fragments with another precursor's overlapping candidate, raw01 q ≤ 0.05):

  | class | share |
  |---|---|
  | human | 17 % |
  | entrapment | 41 % |
  | decoy | 36 % |

  This matches the old library. In about 95 % of cases the claimant is not a stronger winner.
- **DIA-NN 2.2** with its own decoys, raw01, same settings as `diann-ourlib-20261112`, 13:27 wall time:

  | library | q ≤ 0.01 N | FDP_c | q ≤ 0.05 FDP_c | N@FDP_c 1 % / 2 % |
  |---|---|---|---|---|
  | old 100K | 2357 | 0.0180 | 0.0588 | 2150 / 2392 |
  | PAIRLIB | 2361 | 0.0144 | 0.0598 | 2067 / 2522 |

  - DIA-NN's FDP_c at q ≤ 0.01 drops from 0.018 to 0.014 with parents present. That is within count noise (21 vs 17
    entrapments).
  - At q ≤ 0.05 DIA-NN sits at about 1.2× nominal, and our pipeline at about 1.6×.

## Reading
1. The 1.5× entrapment excess is not mainly an artifact of the 100K sampling. With every parent present it persists
   at every q cut, on both raws and in both arms. RESULT_QCAL's "missing parent" explanation accounts for at most
   about a fifth of it.
2. DIA-NN on the same libraries is closer to nominal (about 1.2× at q ≤ 0.05). Part of the remaining excess is
   specific to our scoring.
3. **Leading hypothesis for the next diagnosis, not yet tested:** the trained scorers learn decoy-specific artifacts
   that entrapments lack.
   - The decoys are pseudo-reversed with a swapped C-terminal K/R. That shifts every y ion, and the decoy keeps its
     target's library intensities.
   - Entrapments are FDRBench shuffles with real-peptide-like ladders.
   - A scorer that separates "decoy-like" rather than "false" ranks equally false entrapments above decoys in the
     tail. Off-peak null windows do not show this: null AUC is about 0.497–0.499.
   - Evaluation-only test: among windows below the q ≤ 0.20 cut, how well do the label-free features separate decoy
     winners from entrapment winners? An AUC well above 0.5 means decoy artifacts are present.
4. The interference graph is still worth building for ranking: shared signal marks false winners about 2.3× more
   often than human winners. On its own it will not fix calibration, because decoys and entrapments carry shared
   signal at similar rates.
