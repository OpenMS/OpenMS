# RESULT_MT_H: PLAN_MT on held-out H (20 Nov 2026, reported, no criteria)

H (22474.d) was submitted at 21:18:11Z, after the D verdict was committed in docs `666963d271` at 21:18:09Z. It used the same container, library, evaluator and arms as D.

## Per seed

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Wall (s) | Peak RSS (GB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CUR | 0 | 30,003 | 367 | 1.22 % | 1.40 % | 29,236 | 684 | 26.4 |
| CUR | 1 | 30,119 | 346 | 1.15 % | 1.31 % | 29,286 | 737 | 26.6 |
| CUR | 2 | 30,084 | 356 | 1.18 % | 1.35 % | 29,434 | 684 | 25.9 |
| CUR | 3 | 29,996 | 340 | 1.13 % | 1.29 % | 29,353 | 687 | 25.4 |
| CUR | 4 | 29,586 | 316 | 1.07 % | 1.22 % | 29,209 | 626 | 26.3 |
| MT | 0 | 30,749 | 368 | 1.20 % | 1.37 % | 29,857 | 645 | 27.4 |
| MT | 1 | 30,547 | 349 | 1.14 % | 1.30 % | 29,873 | 645 | 27.5 |
| MT | 2 | 30,706 | 359 | 1.17 % | 1.33 % | 29,925 | 648 | 27.3 |
| MT | 3 | 30,800 | 355 | 1.15 % | 1.32 % | 30,064 | 665 | 27.9 |
| MT | 4 | 30,607 | 338 | 1.10 % | 1.26 % | 29,814 | 711 | 29.2 |

## Pooled

| Arm | Lower (90 %) | Combined (90 %) | Median N @ FDP_c 1 % | Pair win |
|---|---|---|---:|---:|
| CUR | 1.15 % (1.11-1.20) | 1.31 % (1.26-1.37) | 29,286 | 0.501 |
| MT | 1.15 % (1.11-1.20) | 1.32 % (1.27-1.37) | 29,873 | 0.517 |

## Reading

- **The gain is consistent but smaller on H: +2.0 %** (29,286 to 29,873), and every MT seed is above every CUR seed. On D it was +6.5 %. H already reaches 71 % of DIA-NN (41,272 at FDP_c 1 %), so there are relatively fewer faint right peaks left for the block to recover. It is now 72 %.
- **Calibration on H is unchanged.** The pooled lower bound stays at 1.15 % and the combined bound goes from 1.31 % to 1.32 %. The D bound increase (+0.10 pp) does not reproduce on H.
- **Pair fairness moves from 0.501 to 0.517.** That is the same direction as D (0.533), but smaller, and well under the 0.56 plan limit.
- Cost: median wall 648 s against 684 s. Peak RSS is about +1.2 GB (27.5 against 26.3 GB), from the 9 extra features per window.

## Decision

D passed every frozen criterion, and H confirms the direction with no calibration cost. Following the PLAN_MT outcome rule, `-matched_filter` is made the default in OpenMS `01c5891ce0`. On r ≈ 1.56 D seed 0 the new default is byte-identical to the tested `-matched_filter` arm, and `-no_matched_filter` is byte-identical to the previous default (`precursors.tsv`, run table, `calibration.json` and `features.tsv`). 25/25 ctests pass.
