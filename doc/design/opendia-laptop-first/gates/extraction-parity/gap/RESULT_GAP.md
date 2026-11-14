# Where OpenDIA loses DIA-NN's IDs (arm E, 14 Nov 2026)

Diagnostic only. DIA-NN 2.2 (same 6-fragment targets, its own decoys) is used as a reference annotation, not as ground truth. Labels choose the groups. Nothing here is fitted into a method. Scripts and outputs are in this directory. The data are in `OpenMS-r1-benchmark-20261107/gap-20261114/`.

## 1. Funnel (human DIA-NN q <= 0.01 IDs, arm E D11FDWT, N@FDP_c 1 % cut)

| Stage | raw01 | raw02 |
|---|---|---|
| DIA-NN human IDs in our library | 2344 | 2309 |
| Lost: not in the scored population | 24 | 26 |
| Lost: no candidate peak within 6 s of DIA-NN's apex | 75 | 75 |
| Lost: candidate exists but is not our winner | 45 | 36 |
| **Lost: right peak is the winner, below the cut** | **626** | **644** |
| Accepted with DIA-NN's peak | 1574 | 1528 |
| We accept, DIA-NN does not (q <= 0.01) | 114 | 95 |

- **Search and peak choice are nearly closed.** 94 % of DIA-NN IDs have DIA-NN's peak as our winner. The search/RT losses total 99 / 101. Of the 24 population losses on raw01, 79 % lie outside the 90.8 s RT window. Their median map error is 153 s.
- **About 85 % of the gap is ranking.** The right peak is found and chosen, then scored below decoys. This matches the earlier funnels (stage D, V7/V16).

## 2. What the rejected 626 look like (raw01, `REJECTED_*.txt`, `DIAG2_*.txt`)

- **They are weak, not absent.** The median OpenDIA apex intensity is 8.3k against 32.8k for accepted IDs. DIA-NN's precursor quantity is about 7x lower, and its Evidence is 2.62 against 3.85.
- **DIA-NN is confident in them anyway.** The median q is 0.0034, the sum of PEP is 11.0 (about 11 expected false), and only 7 have PEP > 0.05.
- **Acceptance falls with abundance.** Acceptance of right-peak winners by DIA-NN quantity quintile, low to high: 0.40 / 0.58 / 0.76 / 0.89 / 0.94.
- **Fewer fragments, less coherent.** Fragments with >= 3 nonzero scans in the central +-5.6 s (median): 3 against 6. Fragment co-elution: 0.30 against 0.74. 12-fragment support: 1 against 6.
- **They are not marginal.** The median score is 7.2 against a cut of 25.8, and only 2 % lie within 1 unit of the cut. Rescoring has to move them a long way, not just nudge the threshold.

## 3. Is the missing evidence already in our arrays? (`EVIDENCE_BAND_*.txt`)

Each rejected ID is matched to the 3 decoy winners nearest in score. As a specificity control, entrapment winners in the same score band get their own matched decoys.

| Label-free feature (central +-5.6 s) | AUC rejected > decoy | AUC entrapment > decoy (band >= R median, n 128) |
|---|---|---|
| fragment ppm SD across fragments (lower is better) | **0.279 (0.263–0.297)**, so 0.72 inverted | 0.460 (0.415–0.510) |
| top fragment share of area (lower is better) | **0.315**, so 0.69 inverted | 0.506 |
| abs(weighted mean ppm) (lower is better) | **0.339**, so 0.66 inverted | 0.464 |
| MS1 M+0 vs fragment-sum correlation | 0.598 | 0.507 |
| MS1 M+1 vs M+0 correlation | 0.559 | 0.532 |
| fragments with >= 3 nonzero scans | 0.559 | 0.536 |

- **Mass accuracy across fragments is the strongest unused signal.** At matched score, real weak peptides have fragments at a consistent calibrated m/z. Decoy winners are mixtures of unrelated ions at scattered m/z. Entrapment in the same band behaves like decoys (CI includes 0.5), which is the specificity property we need.
- **One-fragment dominance** is the second signal. Decoy winners are often one strong interfering ion plus noise.
- **MS1 co-elution** helps a little. Earlier tests of MS1 features (D11M) raised FDP_c, so it needs care.
- **Caveat.** V5 tested a related ppm feature set (featPC: tight-fragment count, median abs ppm) and it failed rule 2 on raw01 only. Its features were medians per window. The cross-fragment SD and the dominance share are new.

## 4. Diagnostic ceiling (`CEILING*_raw01.txt`)

A labelled HistGradientBoosting model was trained on the 2200 right-peak positives against all decoy winners, with 5-fold CV. This is an upper-bound probe, not a method. The table gives TPR among the rejected 626.

| Feature set | at today's decoy rate | at 4x that rate |
|---|---|---|
| D11FDWT (current) | 0.008 (by definition) | 0.312 |
| 20 scorer inputs | 0.165 | 0.431 |
| + 9 evidence features | 0.217 | 0.505 |
| raw arrays only (6 XIC, 4 MS1, lib) | 0.010 | 0.126 |

- **About 100–140 of the 626 are recoverable at today's decoy rate** with current evidence, even with labels. A labelled model is optimistic, and on all positives it does worse than D11FDWT at the 1x point.
- Even with labels, most of the 626 stay below the decoys. **The signal that separates them from decoys is mostly not in the 6-fragment, +-15 s, 1-cycle-binned arrays.**

## 5. What would recover more (ranked by expected yield; each needs its own frozen plan)

1. **Fragment mass consistency plus dominance features** (section 3).
   - These are label-free, scale-free and computed from arrays we already have (ppm.npy, xic.npy).
   - Expected gain is modest: 50–150 at 1 %.
   - Cheapest next arm. Base D11FDWT on arm E.
2. **More fragments at search time, not just at feature time.**
   - DIA-NN also uses 6 here. But the rejected IDs have only 3 fragments with signal, and the 12-fragment features already exist (f12_*; support 1 vs 6).
   - The real lever is choosing fragments the instrument actually sees. A run-specific fragment re-ranking from anchors (observed vs predicted intensity per ion type and position) is label-free.
3. **Interference removal before scoring** (`DESIGN_INTERFERENCE_GRAPH.md`).
   - Decoys win the matched band by borrowing one strong ion (top_share).
   - Shared-fragment assignment across co-eluting candidates (EPIFANY-style) would take that ion away from the decoy, or from the weaker claimant.
   - V13 showed candidate-ablation SSE separates stage-D from decoys at AUC 0.69, with entrapment at 0.53.
   - This is the largest structural lever left, but it is the most work.
4. **Score calibration with the new calibrations.**
   - The scoring chain has not been rerun with commits 0e4a234762 / 20f9f12335.
   - Free to try, but it changes only RT/ppm windows, and the RT/search losses total about 100.

**Follow-up:** steps 1–3 were run as PLAN_GAP (`arms/RESULT_GAP_ARMS.md`).
- The explaining-away block (I) passes. It recovers 149 / 234 of the rejected IDs, and N@FDP_c 1 % rises to 1879 / 1921.
- A and G fail.

Not worth pursuing: the RT window (about 4 % of losses), peak choice (2 %), and nudging the threshold (the rejected IDs sit 18 score units below the cut).
