# RESULT_LIBSIZE: the library-size penalty is no longer the main gap (21 Nov 2026, diagnostic)

The precheck and decision rule are in `PRECHECK_LIBSIZE.md`, committed (`9c286fe35c`) before any job ran. All 18 jobs (6 DIA-NN, 12 OpenDIA) exited 0. The evaluation is `libsize.py eval` (sha256 `67bc1334…`), and the full output is in `EVAL_LIBSIZE.txt`.

## Libraries built

| Level | Targets | Entrapment | r | Decoys (OpenDIA) | DIA-NN transitions |
|---|---:|---:|---:|---:|---:|
| r0 | 172,608 | 0 | 0 | 172,505 | 1,018,701 |
| r1 | 345,118 | 172,510 | 0.999 | 344,725 | 1,970,213 |
| r3 | 690,626 | 518,018 | 3.001 | 689,713 | 3,878,190 |
| r7 (existing) | | | 7.054 | | |

## HeLa IDs at each tool's own q ≤ 0.01 (OpenDIA is the mean of seeds 0 and 1)

| Run | Level | OpenDIA | DIA-NN | OpenDIA / DIA-NN | OpenDIA share of r0 | DIA-NN share of r0 |
|---|---|---:|---:|---:|---:|---:|
| D | r0 | 19,069 | 27,901 | **0.683** | 1.000 | 1.000 |
| D | r1 | 17,964 | 27,941 | 0.643 | 0.942 | 1.001 |
| D | r3 | 16,676 | 26,434 | 0.631 | 0.875 | 0.947 |
| D | r7 | 15,781 | 25,881 | 0.610 | **0.828** | **0.928** |
| H | r0 | 34,735 | 44,866 | **0.774** | 1.000 | 1.000 |
| H | r1 | 33,226 | 43,463 | 0.764 | 0.957 | 0.969 |
| H | r3 | 32,406 | 42,426 | 0.764 | 0.933 | 0.946 |
| H | r7 | 30,839 | 41,231 | 0.748 | **0.888** | **0.919** |

N at combined FDP 1 %, OpenDIA / DIA-NN: D r1 0.618, r3 0.616, r7 0.593. H r1 0.750, r3 0.754, r7 0.741. OpenDIA's combined bound: D r1 1.11 %, r3 0.99 %, r7 1.25 %. H r1 1.04 %, r3 1.19 %, r7 1.33 %.

## Decision rule

| Rule | D | H | Verdict |
|---|---|---|---|
| 1. Excess penalty X = P_OpenDIA − P_DIA-NN at r7 (threshold: D ≥ 15 pts, H ≥ 10 pts) | 17.2 − 7.2 = **10.0 pts** | 11.2 − 8.1 = **3.1 pts** | **Does not hold** |
| 2. OpenDIA / DIA-NN at r0 (threshold: D ≤ 0.75, H ≤ 0.85) | **0.683** | **0.774** | **Holds** |
| 4. Validity: OpenDIA combined bound at r1 and r3 ≤ 1.5 % | 1.11 %, 0.99 % | 1.04 %, 1.19 % | Valid |

**Verdict: rule 2.** Most of the gap is already there with no entrapment at all. As a share of the log-ratio at r7, the r0 gap accounts for 77 % on D (ln 0.683 / ln 0.610) and 88 % on H (ln 0.774 / ln 0.748). The earlier "−49 % versus −10 %" penalty belonged to older code. On the current default, OpenDIA loses 17 % (D) and 11 % (H) going from r0 to r7, against 7 % and 8 % for DIA-NN. A real but secondary excess remains on D (10 points).

## Overlap and retention (seed 0)

- **r0 overlap:** OpenDIA's HeLa IDs are almost a subset of DIA-NN's (92.0 % on D, 93.7 % on H). DIA-NN finds 37 % (D) and 27.5 % (H) of its IDs that OpenDIA misses, even on the plain library.
- **Retention of r0 IDs at r7:** OpenDIA keeps 0.816 (D) and 0.876 (H). DIA-NN keeps 0.861 and 0.895. The retention curves match up to r3 and split only between r3 and r7 on D.
- **Score threshold climb:** the OpenDIA score threshold at q = 0.01 rises from 6.2 to 13.5 on D (5.3 to 9.8 on H). The decoy tail grows with the number of competing targets, which is the expected mechanism for the remaining D excess.

## What this means for the next steps

- **Step 2 (scorer ceiling)** is judged at r0 first, as rule 2 requires. The r0 gap (DIA-NN finds about 8,800 more HeLa IDs on D and 10,100 more on H) is an extraction or scoring problem at a fixed search space, not an FDR-competition problem.
- The B1/C funnel and the matched-filter and RT gains were all measured at r7. They should be re-read as mostly r0 effects.
- The D-only excess (10 points, mostly between r3 and r7) is a secondary target. It fits a crowded run where entrapment peptides find chance co-eluting signal. It does not justify its own lever until the r0 gap is understood.
