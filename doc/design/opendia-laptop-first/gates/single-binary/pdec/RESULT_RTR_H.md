# RESULT_RTR_H: PLAN_RTR on held-out H (21 Nov 2026, reported, no criteria)

H (22474.d) was submitted at 07:32:10Z, after the D verdict was committed in docs `09285d4462` at 07:32:08Z. It used the same container, library, evaluator and arms as D. The CUR arm's `precursors.tsv` is byte-identical to the PLAN_MT `mt` arm on all 5 H seeds.

## Per seed

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % | Wall (s) | Peak RSS (GB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| CUR | 0 | 30,749 | 368 | 1.20 % | 1.37 % | 29,857 | 682 | 27.4 |
| CUR | 1 | 30,547 | 349 | 1.14 % | 1.30 % | 29,873 | 675 | 27.5 |
| CUR | 2 | 30,706 | 359 | 1.17 % | 1.33 % | 29,925 | 694 | 28.5 |
| CUR | 3 | 30,800 | 355 | 1.15 % | 1.32 % | 30,064 | 751 | 26.6 |
| CUR | 4 | 30,607 | 338 | 1.10 % | 1.26 % | 29,814 | 675 | 28.0 |
| RTR | 0 | 31,260 | 375 | 1.20 % | 1.37 % | 30,543 | 1,236 | 27.3 |
| RTR | 1 | 31,145 | 352 | 1.13 % | 1.29 % | 30,618 | 1,221 | 28.1 |
| RTR | 2 | 31,387 | 374 | 1.19 % | 1.36 % | 30,665 | 1,422 | 28.7 |
| RTR | 3 | 31,490 | 367 | 1.17 % | 1.33 % | 30,868 | 1,288 | 27.5 |
| RTR | 4 | 31,282 | 359 | 1.15 % | 1.31 % | 30,552 | 1,360 | 27.1 |

**RT refit:** the training size was 30,547-30,800 IDs. The p95 RT residual fell from 90.2-91.0 s to 63.9-64.8 s, and the RT window narrowed from ±139.3 s to ±98.5-99.2 s.

## Pooled

| Arm | Lower (90 %) | Combined (90 %) | Median N @ FDP_c 1 % | Pair win |
|---|---|---|---:|---:|
| CUR | 1.15 % (1.11-1.20) | 1.32 % (1.27-1.37) | 29,873 | 0.517 |
| RTR | 1.17 % (1.12-1.21) | 1.33 % (1.28-1.38) | 30,618 | 0.517 |

## Reading

- **H confirms D's direction, more cleanly: +2.5 %** (29,873 to 30,618), and every RTR seed is above every CUR seed (lowest RTR 30,543, highest CUR 30,064). On D, one seed was flat.
- **Calibration is unchanged on H** (+0.02 / +0.01 pp, well inside the noise) and pair win is identical. On D it improved. Neither run shows any calibration cost.
- **The relative gain is not larger on H,** although H's RT window shrank much more (−40 s against −10 s). H's B1 and B2 losses were already smaller, so a tighter window helps less there than the precheck's residual numbers suggested.
- **Cost is 1.89× wall time** for the second pass, with no memory change.
- **Against DIA-NN on H** (41,272 at FDP_c 1 %), this is 74 %, up from 72 %.

## Decision

D passed every frozen criterion and H confirms it with no calibration cost. Following the PLAN_RTR outcome rule, `-rt_refine` is now the default in OpenMS `14c3c54113`, and `-no_rt_refine` keeps the single pass. On r ≈ 1.56 D seed 0, the new default is byte-identical to the tested `-rt_refine` arm, and `-no_rt_refine` is byte-identical to the previous default, across `precursors.tsv`, the run table, `calibration.json` and `features.tsv`. 26/26 tests pass. The 1.8-1.9× wall time is the trade-off, and `-no_rt_refine` keeps the faster single pass for anyone who needs it.
