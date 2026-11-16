# RESULT_RTPROFILE: ID-based cross-fitted recalibration with an RT-dependent window (design runs, verdict)

Written 16 November 2026, after the raw 01/02 runs and **before** any raw 03 run. The plan is `PLAN_RTPROFILE.md`
(sha256 file present, docs ae41250942).

## Setup (as frozen)
- **Code:** OpenMS f86fb0db22.
  - Per-precursor RT windows in StreamSearch.
  - OpenDIAWorkflow `-id_recalibration`.
  - Tests: two search tests. The multi-chunk test fails when the chunk-range widening is removed (checked).
- **Container:** `opendia-rtprofile-f86fb0db22.sif` (sha256 6937ad84…).
- **Codon:** `/nfs/research/juan/singj/opendia/rtwin2-20261116/rtprofile`, 8 threads, about 10 min per run (two passes).
- **Inputs:** self-calibrating, arm-E `library.pqp` + `library12_full.pqp`, seeds 0-4.
- **Arms:** CUR = `<run>.precursors.pass1.tsv`, NEW = `<run>.precursors.tsv` (pass 2).
- **Evaluation:** `eval_run.py`. The tool's q-values match the independent recompute on every file.

## What pass 2 used
- ID anchors: 1671-1772 per run (raw01) and 1671-1814 (raw02). The calibration had 426 / 546 anchors.
- Held-out map error (σ per fold): 21.4-22.8 s on both raws.
- Window half widths, against T = 151 / 124 s:

  | run | narrowest | widest (late gradient end) |
  |---|---|---|
  | raw01 | 46-50 s | 524-622 s |
  | raw02 | 29-38 s | 304-346 s |

## N@FDP_c 1 % (FDP_c at q ≤ 0.01)
| run | arm | s0 | s1 | s2 | s3 | s4 | median |
|---|---|---|---|---|---|---|---|
| raw01 | CUR | 1815 (.0079) | 1804 (.0082) | 1784 (.0047) | 1728 (.0095) | 1827 (.0046) | **1804** |
| raw01 | NEW | 1773 (.0091) | 1874 (.0088) | 1845 (.0068) | 1854 (.0091) | 1719 (**.0104**) | **1845** |
| raw02 | CUR | 1894 (.0088) | 1928 (.0044) | 1935 (.0033) | 1867 (.0046) | 1909 (.0012) | **1909** |
| raw02 | NEW | 1883 (.0077) | 1939 (.0057) | 1937 (.0066) | 1882 (.0087) | 1889 (.0070) | **1889** |

N@FDP_c 0.5 % / 2 %, medians:

| run | arm | 0.5 % | 2 % |
|---|---|---|---|
| raw01 | CUR | 1616 | 1993 |
| raw01 | NEW | 1540 | 2012 |
| raw02 | CUR | 1800 | 2018 |
| raw02 | NEW | 1680 | 2033 |

## Criteria
| criterion | raw01 | raw02 |
|---|---|---|
| 1. NEW median ≥ CUR − 1 % | 1845 vs 1786: pass (+2.3 %) | 1889 vs 1889.9 needed (−1.05 %): **fail** (1 ID short) |
| 2. NEW FDP_c ≤ 0.010 every seed | seed 4: 0.0104, **fail** | max 0.0087: pass |
| 3. NEW ≥ CUR + 2 % on at least one raw | raw01 +2.3 %: yes | |

**Verdict: FAIL.** Criterion 1 fails on raw02 and criterion 2 fails on raw01. Per the plan, `-id_recalibration` is
removed and nothing is retuned.

## Reading
- Recalibrating on the run's own IDs makes the map at the gradient ends data-driven. The window then follows the error
  profile: 3-4× wider at the late end, 2-4× narrower in the middle.
- The gain shows at looser FDR: +19 / +15 median at 2 %, +41 on raw01 at 1 %.
- At strict FDR, however, the pass-2 ranking is less stable:
  - The 0.5 % medians fall by 76 / 120.
  - Raw01 seed 4 collapses at 0.5 % (681) and crosses FDP_c 0.010 at q ≤ 0.01.
  - A likely cause, not tested: the very wide end windows (up to 4× T) admit many more chance peaks for the end
    precursors and their decoys. The scorer's RT-error feature is in seconds, so it cannot tell a 100 s error in a
    600 s window from one in a 50 s window.
- Raw 03 (held out) follows below, reported descriptively.
