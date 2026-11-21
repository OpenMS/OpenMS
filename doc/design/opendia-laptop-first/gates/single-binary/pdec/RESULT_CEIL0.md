# RESULT_CEIL0: the scorer is a small lever; most missed IDs need new evidence or new candidates (21 Nov 2026, diagnostic)

The precheck and decision rule are in `PRECHECK_CEIL0.md`, committed (`edfd981cfe`) before any job ran. The script is `precheck_ceil0.py` (sha256 `19955add…`). All four `-features_out` runs exited 0, and each `precursors.tsv` is byte-identical to its reference run (r0: `libsize-20261121/r0/opendia/{D,H}_s0`, r7: `entlib-rtr-20261121/runs/{D,H}_rtr_s0`). The full output is in `EVAL_CEIL0.txt`. Each feature table has 45 features.

## HeLa IDs at own q ≤ 0.01 (seed 0)

| Arm | D r0 | H r0 | D r7 | H r7 |
|---|---:|---:|---:|---:|
| DIA-NN | 27,901 | 44,866 | 25,881 | 41,231 |
| A tool | 19,097 | 34,702 | 15,859 | 30,888 |
| B label-free GBM | 20,261 (+6.1 %) | 37,178 (+7.1 %) | 15,961 (+0.6 %) | 31,935 (+3.4 %) |
| Bt label-free, hard decoys | 20,338 (+6.5 %) | 36,919 (+6.4 %) | 16,785 (+5.8 %) | 32,917 (+6.6 %) |
| **C ceiling (DIA-NN labels)** | **20,644 (+8.1 %)** | **37,444 (+7.9 %)** | 15,583 (−1.7 %) | 32,462 (+5.1 %) |
| Ct ceiling, hard decoys | 16,354 (−14.4 %) | 34,957 (+0.7 %) | 15,274 (−3.7 %) | 31,600 (+2.3 %) |
| Tool / DIA-NN | 0.684 | 0.773 | 0.613 | 0.749 |
| Best ceiling / DIA-NN | 0.740 | 0.835 | 0.602 | 0.787 |

At r7, N at combined FDP 1 % and the lower / combined bounds were:

| Arm | D | H |
|---|---|---|
| A | 15,465 (1.12 / 1.27 %) | 30,543 (1.20 / 1.37 %) |
| B | 15,387 (1.21 / 1.39 %) | 31,250 (1.26 / 1.43 %) |
| Bt | **16,131 (+4.3 %)** (1.23 / 1.40 %) | **32,353 (+5.9 %)** (1.22 / 1.39 %) |
| C | 15,330 (1.01 / 1.15 %) | 31,666 (1.24 / 1.41 %) |

## Decision rule

| Rule | D | H | Verdict |
|---|---|---|---|
| G = max(C, Ct) / A − 1 at r0. Rule 1 needs ≥ 15 % (D) or ≥ 10 % (H); rule 2 needs < 5 % on both | **+8.1 %** | **+7.9 %** | **Middle band (rule 3): the funnel decides** |
| 4. B or Bt ≥ +5 % at r0 on both runs, ≥ +3 % at N@FDP_c 1 % at r7 on both, lower bound +0.2 pp or less | Bt: +6.5 %, +4.3 %, +0.11 pp | Bt: +6.4 %, +5.9 %, +0.02 pp | **Bt qualifies** (B fails at r7 on D) |

**Reading.** Even with DIA-NN's own labels, a learner on our 45 features closes only about 8 % of the 46 % (D) and 29 % (H) shortfall at r0, reaching 74 % and 84 % of DIA-NN. The features lack the information for most of DIA-NN's extra IDs. A label-free tail-focused learner (Bt) already gets most of what the ceiling gets. At r7 it beats the ceiling, because DIA-NN-labelled positives include faint windows that sit inside the decoy tail.

## Funnel of DIA-NN-only HeLa IDs (tool score)

| Class | D r0 (10,344) | H r0 (12,360) | D r7 (10,865) | H r7 (11,973) |
|---|---:|---:|---:|---:|
| A_low: the right window is the best, but scores too low | **61.5 %** (6,362) | 40.3 % (4,981) | 62.3 % | 50.9 % |
| B2: the right window loses to another window | 11.4 % (1,179) | 11.9 % (1,471) | 10.1 % | 11.0 % |
| C: no window within 6 s of DIA-NN's RT | 27.2 % (2,814) | **47.8 %** (5,908) | 27.6 % | 38.1 % |

- The class shares are almost flat across DIA-NN quantity quartiles. Even at Q4 (the brightest quarter), C is 30 % on D and 54 % on H. So C is not a faint-signal effect. It points to candidate generation: the search window, the valley merge or the peak picker.
- **TPR of the right window at decoy FPR 1e-4** (tool / ceiling C): D r0 0.395 / 0.373, H r0 0.576 / 0.548. Supervised training raises TPR at FPR 1e-2 (D 0.674 → 0.734) but not at the 1e-4 tail that decides 1 % FDR. This is the same pattern as PRECHECK_SCORER36. The A_low windows look like decoys on every feature we have.
- These class definitions differ from the earlier B1/C funnel. Here C means no window within 6 s of DIA-NN's apex, so the numbers are not directly comparable.

## What follows

1. **Candidate recovery (C) is the largest single class on H** (5,908 at r0) and the second largest on D (2,814). It does not depend on signal strength. This is step 3. The first precheck: for C IDs, find where the right RT falls relative to our search window and the candidates we did emit (outside the window, merged into a neighbour, or below the picking threshold). It needs no new runs, because the feature tables and libraries are already on codon.
2. **New evidence for A_low** (6,362 on D, 4,981 on H). The scorer cannot fix these with current features. They need information we do not compute.
3. **Bt (label-free, hard-decoy training)** qualifies under rule 4 for an r ≈ 1.56 precheck and then a frozen 5-seed plan. It is worth about +4-6 % at r7. Caveat: the same idea in Python (`precheck_scorer36_tailfit.py`, PRECHECK_NOCAND section 3, older default) gave only +1.2 % with a worse lower bound at r ≈ 1.56. So the r ≈ 1.56 check of Bt on the current default comes before any C++ work (PRECHECK_BT156).

No code or default changes.
