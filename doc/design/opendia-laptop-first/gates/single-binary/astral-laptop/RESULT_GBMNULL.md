# RESULT_GBMNULL: the null-region flip test cannot tell when the GBM is unsafe (23 Nov 2026, FAIL)

The plan is `PRECHECK_GBMNULL.md`. It was frozen and pushed in docs `88b994923c` before any statistic or raw 03 run.

- **Fix to the script after freezing (no change to the statistic):** the first pass crashed on a pandas name clash (column `gt` shadowed `DataFrame.gt`). The columns were renamed. The region, test and level are unchanged.
- **Runs:** the Astral GBM runs are the post-fix set `astral-gbmfix-20261123` (the RESULT_ASTRAL_LAPTOP numbers). Three batch-1 directories had crashed.
- **Anchor gate on raw 03: not met.** A0 (`gbmdef`, no flags) crashed on raw03_s0 with `pooledPairedTDC: a fold has zero decoy spread` after `trees 127/0/120` in GBM round 3. That is the bug fixed in `88e5da3023`, which the old `gbmdef` container lacks. By the frozen rule, raw 03 does not count. Its numbers below are informational only.
- **The verdict does not depend on raw 03.** raw 02 alone fails.

## Statistic (label-free; p < 0.01 means skip)

| Set | Decisions (seeds 0-4) | p values |
|---|---|---|
| D | use, use, use, use, **skip** | 0.20, 0.028, 0.026, 0.98, 0.0003 |
| H | use, use, use, **skip**, **skip** | 0.034, 1.0, 0.80, 9e-5, 2e-4 |
| A1 | use ×5 | 0.92, 0.40, 1.0, 1.0, 1.0 (s4: the GBM changed nothing) |
| A2 | use ×5 | 0.18, 1.0, 0.026, 0.63, 1.0 |
| A3 (informational) | use ×5 | 0.51, 0.081, 0.11, 0.40, 0.63 |

In the null region the target win share stays near 0.50 under both scores on every set (0.45-0.505). It never rises on Astral.

## Grading (pooled lower / combined bound at own q ≤ 0.01, median N @ FDP_c 1 %)

| Set | Never (no GBM) | Always (GBM) | Rule | Oracle | |
|---|---|---|---|---|---|
| D | 1.10 / 1.25 %, 19,704 | 1.11 / 1.27 %, 21,396 (admissible) | 1.11 / 1.27 %, 21,054 | 21,396 | pass |
| H | 1.15 / 1.31 %, 34,871 | 1.29 / 1.47 %, 37,342 (admissible) | 1.24 / 1.41 %, 37,339 | 37,342 | pass |
| A1 | 0.48 / 0.96 %, 1,943 | 0.61 / 1.22 %, 2,136 (admissible) | = always | 2,136 | pass |
| **A2** | 0.26 / 0.52 %, 2,285 | 0.56 / 1.11 %, 2,146 (**not admissible**) | = always | 2,285 | **fail** (bounds and N) |
| A3 (informational) | 0.38 / 0.76 %, 2,231 | 0.60 / 1.20 %, 2,246 (not admissible) | = always | 2,231 | (would fail on bounds) |

**Verdict: FAIL.** The rule always used the GBM on Astral, including raw 02 where the GBM is inadmissible. On HeLa it skipped 3 of 10 seeds where the GBM was safe and helpful.

## Reading

- **The null region carries no signal.** When the GBM makes q optimistic on Astral, it does not tilt target-against-decoy wins among pairs that are both clearly false. The extra false targets must enter near or above the threshold, where true targets also sit. As the precheck feared, that region cannot be checked without labels.
- **The HeLa "skips" are noise at this level.** With about 170,000 flips, p < 0.01 catches shifts of about 0.4 %. Three of the 10 HeLa seeds do so in a direction that does no harm.
- **Nothing is retuned.** No other region, level or statistic is tried on these runs. The GBM stays as is, and the Astral looseness remains a known limit.
