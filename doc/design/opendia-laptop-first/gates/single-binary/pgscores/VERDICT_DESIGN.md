# PLAN_PGSCORES design verdict: FAIL (criterion 2)

Written 18 November 2026, after the D design runs and before any held-out run (H or Astral raw 03).

**Correction (10:07, from file timestamps):** The Astral design runs (raws 01/02, NEW seeds 0-4) had in fact finished on
codon by 09:09, before this verdict. They had not been fetched or evaluated, and I wrongly wrote that they were still
queued. They cannot change the verdict, and their result is in RESULT_PGSCORES.md. The held-out timing holds: the raw 03
jobs were submitted at 09:14 and H started at 09:14, both after this commit (09:13:44).

| D (design) | N median (min-max) | agreement with DIA-NN, median |
|---|---|---|
| CUR (tip 963543c7f0) | 16,116 (15,957-16,177) | 0.9299 |
| NEW (`-pg_scores`) | 18,792 (18,585-18,885) | 0.9104 |

- **Criterion 1** (N ≥ 16,600): **pass**, +16.6 %. Every NEW seed is above every CUR seed.
- **Criterion 2** (agreement ≥ 0.9249): **fail**, 0.9104.
- **Seed 0 detail** (`AGREE_D.txt`):
  - CUR: 14,910 of 16,021 are in DIA-NN's 1 % set (1,111 are not).
  - NEW: 16,979 of 18,639 (1,660 are not).
  - Of the 2,618 added IDs, 2,069 are DIA-NN IDs and 549 are not, a marginal agreement of 0.79.
  - For context, OpenSWATH SVM on the same run has agreement 0.886 (19,132 IDs, 2,179 not in DIA-NN).
- **Flag off:** D seed 0 output is byte-identical to the tip.
- **Cost:** extraction 56 → 117 s on D, peak memory unchanged (4.4 GB).

**Outcome by the plan:** fail. The code is to be removed and nothing is retuned. The held-out H and raw 03 runs, and the
Astral design runs, are still reported.

**Caveat:** diaPASEF has no entrapment, so this criterion is a proxy. A drop in agreement can mean false IDs, or true IDs
that DIA-NN also misses. The plan fixed the proxy in advance, and the result is reported against it as written.
