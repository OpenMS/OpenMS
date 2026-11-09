# Single TOPAZ scorer: result (partial)

Rule: `PLAN.md` (sha256 7737c2a7, written 2026-11-09T17:05:50Z).

## Default equivalence

The HeLa `crossfit.tsv` is IDENTICAL with the new binary (`hela-default-cmp/CMP.txt`).

## raw 01 S0: FAIL (did not complete, which breaks rule 4)

| Stage | Target precursors at q ≤ 0.01 (decoy-count estimate) |
|---|---|
| Score iteration 0 (seed training) | 347 |
| Score iteration 1 | 149 |
| Iteration 2 | crashed: "no positive seeds outside held-out fold 2" (all 149 positives are in one fold) |
| Base: CNN + LDA, native byion | 1236 at stage-2 q ≤ 0.01; N@1% 1029 |

Run cost: 16:07 wall time, 25.5 GB peak RSS. There is no `crossfit.tsv`, so N@1% and FDP_c cannot be read.

- **Validation loss:** per fold, 0.0037 / 0.0002 / 0.0044 in round 0, falling to 0.00002 to 0.0002 in round 1. Held-out seed AUC is 1.0. With the heuristic features, the seed task is trivially separable, so the net saturates and the held-out decoy tail decides q.
- **Semi-supervised rounds:** they made it worse (347 → 149), the same direction as earlier ss runs.
- **Chain:** it stopped on the failure. S1 and raw 02 did not run.

## Round 2 (`PLAN2.md`, sha256 9073a344): A1 FAIL on both raws (`EVALUATION2.txt`)

There were no score iterations, so q comes from the round-0 TOPAZ logit.

| | raw 01 A1 (augmented) | raw 01 C0 (control) | raw 02 A1 | base raw 01 / 02 (CNN + LDA) |
|---|---|---|---|---|
| N at FDP_c 0.5 / 1 / 1.5 / 2% | 426 / **488** / 552 / 609 | 13 / **13** / 13 / 110 | 434 / **449** / 684 / 736 | 1029 / 1298 at 1% |
| q ≤ 0.01: N, FDP_c | **0** | 347, 0.012 | **0** | 1236 at stage-2 q (raw 02) |
| q ≤ 0.05: N, FDP_c | 425, 0.043 | 830, 0.044 | 566, 0.053 | FDP_c 0.065 / 0.074 |
| R at decoy median, per fold | 0.92 / 0.91 / 0.92 | 0.93 / 0.92 / 0.92 | 0.92 / 0.91 / 0.92 | |
| Held-out seed AUC | 1.0 | 1.0 | 1.0 | |
| Time, peak RSS | 10.8 min, 24.1 GB | 12.7 min, 24.1 GB | 13.5 min, 23.9 GB | |

**Rule check:** rule 1 (N@1%) fails for A1 on both raws. Rules 2 and 3 pass trivially for A1 (zero IDs at q ≤ 0.01), and rule 4 passes.

**Augmentation effect (reported, not a pass condition):** A1's N@1% is 488 against C0's 13. That meets the predeclared ≥ 1.05 × C0 criterion, but the comparison is not clean.
- C0 ranks a few entrapment precursors at the very top, so its FDP_c curve breaks after 13. A1 removed them.
- In exchange, A1's decoy tail got worse: q ≤ 0.01 went 347 → 0 and q ≤ 0.05 went 830 → 425.
- The mixed windows (s × target + decoy) are labelled positive, so the head learns that decoy-like window content can be positive. That plausibly raises the top decoy scores.
- So augmentation helps at the very top but not overall. A stronger version of the same augmentation is not expected to fix the decoy tail.

## Round 3 (`PLAN3.md`, sha256 18f58e92): N1 decoy negatives FAIL on both raws (`EVALUATION3.txt`)

HeLa default `crossfit.tsv` IDENTICAL (`hela-default-cmp3/`).

| | raw 01 N1 | raw 01 C0 (no negatives) | raw 02 N1 | base raw 01 / 02 |
|---|---|---|---|---|
| N at FDP_c 0.5 / 1 / 1.5 / 2% | 739 / **950** / 1003 / 1047 | 13 / 13 / 13 / 110 | 784 / **893** / 1002 / 1085 | 1029 / 1298 at 1% |
| q ≤ 0.01: N, FDP_c | 340, 0.000 | 347, 0.012 | 859, **0.016** | |
| q ≤ 0.05: N, FDP_c | 1420, **0.091** | 830, 0.044 | 1460, 0.073 | FDP_c 0.065 / 0.074 |
| R per fold | 0.95 / 0.93 / 0.95 | 0.93 / 0.92 / 0.92 | 0.92 / 0.95 / 0.93 | |
| Time, peak RSS | 18.6 min, 24.1 GB | 12.7 min | 21.5 min, 23.4 GB | |

**Rule check:**
- Rule 1 (N@1%) fails on both raws.
- Rule 3 fails on raw 01 (0.091 > 0.070).
- Rule 2 fails on raw 02 (0.016 > 0.015).
- Rule 4 passes.

**Predeclared effect criterion (q ≤ 0.01 ≥ 694 on raw 01): not met** (340 against C0's 347). By that criterion, the hypothesis is not supported.

The ranking still changed a lot. N@1% went from 13 to 950 on raw 01, and to 893 on raw 02. That is 92% / 69% of the LDA base, far above every earlier single-scorer arm. The decoy negatives fixed the target ranking (entrapment at the top) but not the q calibration:
- Raw 01 q ≤ 0.05 is anti-conservative (FDP_c 0.091).
- At q ≤ 0.01, raw 01 is overly conservative (340 with FDP_c 0.000, against 950 truly at FDP_c 1%).
- The decoy winners and false targets in the upper tail are still scored differently.
- Held-out seed AUC is still 1.0.

The single scorer is closer, but it is not yet a replacement. CNN + LDA stays the native design.

## Round 4 (`PLAN4.md`, sha256 8d164e55): L1 linear head FAIL on both raws (`EVALUATION4.txt`)

HeLa default `crossfit.tsv` IDENTICAL (`hela-default-cmp4/`).

| | raw 01 L1 | raw 01 N1 | raw 02 L1 | raw 02 N1 | base raw 01 / 02 |
|---|---|---|---|---|---|
| N at FDP_c 0.5 / 1 / 1.5 / 2% | 498 / **553** / 660 / 692 | 739 / 950 / 1003 / 1047 | 485 / **495** / 564 / 668 | 784 / 893 / 1002 / 1085 | 1029 / 1298 at 1% |
| q ≤ 0.01: N, FDP_c | 235, 0.000 | 340, 0.000 | 511, **0.020** | 859, 0.016 | |
| q ≤ 0.05: N, FDP_c | 909, 0.064 | 1420, 0.091 | 861, 0.037 | 1460, 0.073 | 0.065 / 0.074 |
| Held-out seed AUC | 1.0 | 1.0 | 1.0 | 1.0 | |
| Time | 19.8 min | 18.6 min | 17.4 min | 21.5 min | |

**Rule check:**
- Rule 1 fails on both raws.
- Rule 2 fails on raw 02 (0.020).
- Rules 3 and 4 pass.

**Predeclared calibration criterion:** |FDP_c at q ≤ 0.05 − 0.05| shrank on both raws: 0.041 → 0.014 on raw 01, 0.023 → 0.013 on raw 02. So it is **met**.

The cost is ranking power. N@1% fell by 42% (raw 01) and 45% (raw 02), and raw 02 q ≤ 0.01 is still anti-conservative.

The trade-off between ranking (MLP head) and calibration (linear head) is inside one TOPAZ logit. Neither head gets both. The held-out seed AUC of 1.0 in every arm means the seed labels never limit the fit, so the limit is in the score scale that the training objective produces, not in the head width.

## Fold diagnostic (`fold_diag.py`, reported only, no rule)

Pooled TDC standardises each fold by its decoy-winner median and MAD.

Target winners at pooled q ≤ 0.01, per fold:

| Arm | Raw 01 | Raw 02 |
|---|---|---|
| N1 | 296 / 1 / 43 | 208 / 359 / 292 |
| L1 | 107 / 128 / 0 | 252 / 217 / 42 |

The folds are badly unbalanced: raw 01 N1 has 1 in fold 1. Median/MAD scaling does not align the per-fold upper tails of the TOPAZ logit.

Fold-wise TDC on the raw logit does not fix this either:

| Arm | q ≤ 0.01 (N, FDP_c) raw 01 | raw 02 | q ≤ 0.05 FDP_c raw 01 | raw 02 |
|---|---|---|---|---|
| N1 | 822, 0.022 | 652, 0.019 | 0.090 | 0.076 |
| L1 | 337, 0.018 | 331, 0.024 | 0.064 | 0.053 |

So there are two problems:
- between folds, the score scales differ;
- within each fold, the false targets outscore the decoy tail. This is anti-conservative at 1% for every arm.

A per-fold rescaling alone does not fix it.

## Decoy-ness diagnostic (`decoyness.py`, `lda_on_n1.py`; reported only)

This measures whether a score separates two kinds of false winners: entrapment targets and decoys (in-window precursor winners, per fold). TDC is only valid when the AUC is about 0.5.

| Scorer | raw 01 AUC per fold | raw 02 AUC per fold |
|---|---|---|
| Stage-2 LDA (native base) | 0.501 / 0.498 / 0.500 | 0.501 / 0.500 / 0.500 |
| TOPAZ C0 | 0.519 / 0.526 / 0.525 | |
| TOPAZ N1 | 0.518 / 0.521 / 0.518 | 0.525 / 0.515 / 0.524 |
| TOPAZ L1 | 0.524 / 0.524 / 0.523 | 0.519 / 0.524 / 0.525 |
| TOPAZ A1 | 0.523 / 0.525 / 0.525 | 0.523 / 0.527 / 0.528 |

**Every TOPAZ logit scores false targets above decoys** (AUC 0.52). It has learned decoy-specific features. This is consistent with its anti-conservative q at 1%, but it is not a full explanation: the LDA has a global AUC of 0.50 and is still anti-conservative at 5% (FDP_c 0.065 / 0.074). The global AUC does not resolve the upper tail.

The LDA is at 0.50, even though one of its inputs is the CNN logit. A linear refit on decoy winners is what removes the decoy-ness.

Offline check: a fold-crossed 3-round LDA on the N1 logit plus the 11 native heuristics gives AUC 0.50 again:

| Raw | N@1% | q ≤ 0.01 (N, FDP_c) | q ≤ 0.05 (N, FDP_c) |
|---|---|---|---|
| raw 01 | 959 (base 1029) | 1106, 0.018 | 1720, 0.091 |
| raw 02 | 1218 (base 1298) | 1200, 0.010 | 1709, 0.066 |

- It is below the grid2d CNN + LDA base on both raws.
- The q ≤ 0.05 miscalibration of raw 01 stays. It is the known stage-2 excess (CALIB §1–12), not TOPAZ-specific.

**Reading:** the TOPAZ trace network learns a decoy artefact that a linear stage over decoy winners removes. With that stage, TOPAZ N1 is still not better than grid2d. The remaining single-scorer route has to stop the network from learning decoy-ness, for example by making the decoy generator match more target features (§12) or by an adversarial or decoy-matched training objective. That is a design choice for the user.

## Decision (rounds 1–4)

CNN + LDA stays the native design. `-topaz_features`, `-augment_*`, `-topaz_decoy_negatives` and `-topaz_linear_head` stay opt-in and experimental.
- Every arm fits its about 600 seeds per fold perfectly (held-out seed AUC 1.0).
- Decoy negatives (N1) fixed the ranking. The linear head (L1) fixed the calibration at 5% but lost about 40% of N. No single TOPAZ logit gave both.
- The LDA reaches both by refitting a linear score on about 177k decoy winners per fold, with semi-supervised target positives. Its features include the CNN logit.
- Next options:
  - more training signal, for example MIL over all precursors with the heuristic features (bounded by the arm D budget);
  - a strictly regularised head, for example a linear head on the features plus the TOPAZ embedding, to keep the LDA's low capacity.
