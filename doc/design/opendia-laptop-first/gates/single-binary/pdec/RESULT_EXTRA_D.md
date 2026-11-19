# RESULT_EXTRA_D: PLAN_EXTRA on D (19 Nov 2026)

**The plan.** PLAN_EXTRA was frozen in 5df717133a at 18:10:44Z (`PLAN_EXTRA.sha256`, `PLAN_EXTRA.written_at`).
- The D runs were submitted at 18:14:33Z: 10 jobs, 36459118-36459139.
- Container `opendia-pd8-83a0f6c14a.sif` (sha256 c6bcac1b…).
- Library: the PLAN_ODP predicted library (sha256 76f2a6c2…).
- Evaluator: `eval_pdec.py` (sha256 bb03387d…). Output is in `EVAL_EXTRA_D.txt`.

**The baseline reproduces.** CUR seed 0 `precursors.tsv` is byte-identical to PLAN_ODP CUR s0. The pooled CUR numbers
match RESULT_ODP_D exactly (E 458, N 43,243, median N@FDP_c 1 % 8,233).

| | CUR (default) | NEW (`-extra_scores`) |
|---|---|---|
| N at own q ≤ 0.01, per seed | 8,484 / 8,851 / 8,525 / 8,653 / 8,730 | 11,905 / 11,993 / 11,868 / 11,986 / 11,823 |
| **median N at FDP_c 1 %** (range) | **8,233** (7,968-8,286) | **11,332** (11,258-11,542) |
| pooled lower bound at q ≤ 0.01 (90 %) | 0.0106 (0.0098-0.0114) | 0.0116 (0.0109-0.0124) |
| pooled combined bound at q ≤ 0.01 (90 %) | 0.0121 (0.0112-0.0131) | 0.0133 (0.0125-0.0141) |
| P(entrapment beats own decoy, best q ≤ 0.05) | 0.534 | 0.515 |
| wall / peak RSS per run (codon, 8 threads, cached library) | 8:35-11:13 / 37.3-37.6 GB | 8:50-9:20 / 37.6-38.2 GB |

## Criteria
| # | criterion | limit | NEW | verdict |
|---|---|---|---|---|
| 1 | median N at FDP_c 1 % ≥ 1.10 × CUR | ≥ 9,056.3 | **11,332** (1.376×) | pass |
| 2a | lower bound ≤ CUR + 0.0025 | ≤ 0.013091 | **0.011632** | pass |
| 2b | combined bound ≤ CUR + 0.0030 | ≤ 0.015093 | **0.013281** | pass |
| 3 | pair win ≤ 0.56 | ≤ 0.56 | **0.515** | pass |

**Verdict: PASS (all criteria).** Every NEW seed is above every CUR seed at equal true FDR. The worst NEW seed has
11,258 IDs and the best CUR seed 8,286.

## What it means
- **IDs at 1 % true FDP on D go from 8,233 to 11,332 (+37.6 %).** That is 32 % → 44 % of DIA-NN's 26,006.
- **Calibration loosens a little but stays inside the frozen tolerance.**
  - Own q ≤ 0.01 has a true FDP of 1.16-1.33 % (NEW), against 1.06-1.21 % (CUR).
  - Both still fail the absolute PLAN_ODP limits (lower ≤ 1.0 %, combined ≤ 1.2 %), as the plan stated.
- **Pair fairness improves** (0.534 → 0.515). The new evidence does not favour false targets over their decoys.
- **Cost is negligible:** about +0 to 5 % wall time and +0.5 GB memory.

## Next
Per the plan, H (22474.d, seeds 0-4, both arms) is submitted after this verdict commit and reported either way. Making
`-extra_scores` the default is a separate commit with ctests, after H.
