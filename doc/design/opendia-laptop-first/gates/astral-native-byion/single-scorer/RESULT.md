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
