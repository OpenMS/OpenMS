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

## Decision

The single TOPAZ scorer fails on raws 01 and 02, with and without augmentation and with and without score iterations. CNN + LDA stays the native design. `-topaz_features` and `-augment_*` stay opt-in and experimental.
- The shared failure is a saturated fit on about 600 seeds per fold (held-out seed AUC 1.0 in every arm) plus a heavy held-out decoy tail.
- Next options:
  - more training signal, for example MIL over all precursors with the heuristic features (bounded by the arm D budget);
  - a strictly regularised head, for example a linear head on the features plus the TOPAZ embedding, to keep the LDA's low capacity.
