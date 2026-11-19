# RESULT_IM2D_D: PLAN_IM2D on D (19 Nov 2026)

**Plan.** PLAN_IM2D was frozen in 38273668e7 at 20:50:11Z (`PLAN_IM2D.sha256`, `PLAN_IM2D.written_at`).
- The D jobs (36738192-36738203) were submitted at 20:50:25Z.
- Container `opendia-pd9-07131138b8.sif` (sha256 312a8eb3…).
- Evaluator `eval_pdec.py` (sha256 bb03387d…). Output is in `EVAL_IM2D_D.txt`.

**Baseline reproduces.** CUR seed 0 `precursors.tsv` is byte-identical to PLAN_EXTRA NEW seed 0. The CUR pooled
numbers match RESULT_EXTRA_D NEW exactly (E 693, N 59,575, median 11,332).

| | CUR (default) | NEW (`-im_peaks 8`) |
|---|---|---|
| search candidates | 10.44 M | 14.77 M (+4.33 M band peaks, +41 %) |
| N at own q ≤ 0.01, per seed | 11,905 / 11,993 / 11,868 / 11,986 / 11,823 | 11,826 / 12,343 / 11,916 / 12,173 / 12,234 |
| **median N at FDP_c 1 %** (range) | **11,332** (11,258-11,542) | **11,420** (11,173-11,578) |
| pooled lower bound at q ≤ 0.01 (90 %) | 0.0116 (0.0109-0.0124) | 0.0119 (0.0112-0.0127) |
| pooled combined bound at q ≤ 0.01 (90 %) | 0.0133 (0.0125-0.0141) | 0.0136 (0.0128-0.0145) |
| pair win | 0.515 | 0.518 |
| median wall / max peak RSS | 9:09 / 36.8 GB | 13:38 / 49.5 GB |

## Criteria
| # | criterion | limit | NEW | verdict |
|---|---|---|---|---|
| 1 | median N at FDP_c 1 % ≥ 1.03 × CUR | ≥ 11,671.96 | **11,420** (1.008×) | **fail** |
| 2a | lower bound ≤ CUR + 0.0020 | ≤ 0.013632 | 0.011919 | pass |
| 2b | combined bound ≤ CUR + 0.0025 | ≤ 0.015781 | 0.013609 | pass |
| 3 | pair win ≤ 0.56 | ≤ 0.56 | 0.518 | pass |
| 4 | wall ≤ 1.5 × CUR, RSS ≤ 64 GB | ≤ 1.5× / 64 GB | 1.489× / 49.5 GB | pass (wall just inside) |

**Verdict: FAIL (criterion 1).** `-im_peaks` stays opt-in and nothing is retuned.

## What it means
- **The gain is within seed noise.** The NEW and CUR seed ranges overlap (11,173-11,578 against 11,258-11,542).
  Calibration and pair fairness are unchanged.
- **It is expensive on the large library.** The extra candidates (+41 %) cost +49 % wall and +13 GB.
- **Why it does not pay off.** On r ≈ 1.56 (PRECHECK_IM2D) it recovered right-RT candidates for about 9 % of DIA-NN
  IDs, but only 13 % of those were accepted.
  - The bottleneck is scoring faint peaks, not finding them. That holds even more at r ≈ 7, where the extra decoy
    candidates also raise the competition.
  - The +3.5 % seen at r ≈ 1.56 on one seed did not carry over.
- **What would make 2-D picking worth it.** The recovered peaks need evidence that separates them from decoys. A
  mobility-resolved score (features computed on the picked band, not the pooled window) is the natural follow-up. It
  is a separate precheck.

## Next
Per the plan, H is submitted after this verdict commit and reported either way.
