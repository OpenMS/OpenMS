# RESULT_RTWIN2: cross-fitted RT width selection (design runs, verdict)

Written 16 November 2026, after the raw 01/02 runs and **before** any raw 03 run (plan: `PLAN_RTWIN2.md`, sha256 OK).

## Setup (as frozen)
- **Code:** OpenMS `feature/opendia-clean-rewrite` d191f44f08.
  - New code: `CrossFitScorer::selectWidth` / `finalise`, and the OpenDIAWorkflow option `-rt_width_selection`.
  - Unit test: `OpenDIAScorer_test` section 7.
- **Container:** `opendia-rtwin-d191f44f08.sif` (sha256 78ea916e…).
- **Codon:** `/nfs/research/juan/singj/opendia/rtwin2-20261116`, 8 threads, about 6 min per run.
- **Inputs:** self-calibrating. Search library arm E `library.pqp` (md5 e0461beb); scoring library
  `library12_full.pqp` (sha256 218997c0…). Seeds 0-4.
- **Arms:** CUR = `<run>.precursors.full_window.tsv`, RTWIN = `<run>.precursors.tsv`, from the same run.
- **Evaluation:** `gap-20261114/scripts/eval_run.py`; entrapment labels are used only here. The tool's own q-values match
  the independent recompute on every file.

## Widths chosen (half width in s, folds 0/1/2; T = calibrated tolerance)
| run | T | s0 | s1 | s2 | s3 | s4 |
|---|---|---|---|---|---|---|
| raw01 | 150.9 | T / T / 83.0 | 113.2 / T / 83.0 | T / T / 113.2 | T / T / 83.0 | T / T / T |
| raw02 | 123.6 | 68.0 / T / 111.3 | T / 68.0 / T | 74.2 / T / 68.0 | 68.0 / 61.8 / 68.0 | 68.0 / T / 68.0 |

The width is chosen inconsistently across folds and seeds: from 0.55 T to T, for the same run.

## N@FDP_c 1 % (and FDP_c at q ≤ 0.01)
| run | arm | s0 | s1 | s2 | s3 | s4 | median |
|---|---|---|---|---|---|---|---|
| raw01 | CUR | 1815 (.0079) | 1804 (.0082) | 1784 (.0047) | 1728 (.0095) | 1827 (.0046) | **1804** |
| raw01 | RTWIN | 1765 (.0091) | 1768 (.0073) | 1754 (.0048) | 1738 (.0072) | 1827 (.0046) | **1765** |
| raw02 | CUR | 1894 (.0088) | 1928 (.0044) | 1935 (.0033) | 1867 (.0046) | 1909 (.0012) | **1909** |
| raw02 | RTWIN | 1915 (.0088) | 1923 (.0056) | 1898 (.0033) | 1875 (.0078) | 1913 (.0065) | **1913** |

## Criteria
| criterion | raw01 | raw02 |
|---|---|---|
| 1. RTWIN median ≥ CUR − 1 % | 1765 vs 1786 needed (−2.2 %): **fail** | 1913 vs 1890 (+0.2 %): pass |
| 2. RTWIN FDP_c ≤ 0.010 every seed | max 0.0091: pass | max 0.0088: pass |
| 3. RTWIN ≥ CUR + 2 % on at least one raw | no | no |

**Verdict: FAIL.** `-rt_width_selection` stays an option, off by default. Nothing is retuned.

## Reading
- Choosing among 13 widths on two folds' TDC counts is noisy: the choice moves by up to 0.45 T between seeds.
- At the 1 % level, the count it maximises is dominated by seed noise (CUR spans 1728–1827 on raw01). The per-fold
  maximum is therefore a selection on noise, which costs ground truth on raw01 and is neutral on raw02.
- Only the q-value count rises (raw02 s4: 1671 → 1849 at q ≤ 0.01). The true count at FDP_c 1 % barely moves
  (1909 → 1913), and FDP_c rises (0.0012 → 0.0065).
  - So the gain is mostly decoy-tail noise, not better identification.
- Side observation, not part of the criteria: with the full 12-fragment scoring library, self-calibrated CUR raw01 has
  FDP_c ≤ 0.0095 on all five seeds and q ≤ 0.01 counts 1686–1772.
  - Earlier, with `library12_pop`, seed 0 was 1718 at FDP_c 0.0104.

Raw 03 (held out) follows below, reported descriptively against criteria 1-2.

## Raw 03 (held out, run once after the verdict above was committed, docs 8230e2aa3b)
| arm | s0 | s1 | s2 | s3 | s4 | median |
|---|---|---|---|---|---|---|
| CUR | 1907 (.0096) | 1831 (.0058) | 1864 (.0068) | 1913 (.0084) | 1957 (.0076) | **1907** |
| RTWIN | 1863 (.0105) | 1834 (.0057) | 1819 (.0069) | 1915 (.0084) | 1923 (.0087) | **1863** |

Widths: T = 107.4 s; the chosen widths range from 0.55 T to T, varying across folds and seeds.
- Criterion 1: 1863 against the 1888 needed (−2.3 %), **fail**.
- Criterion 2: seed 0 has FDP_c 0.0105, **fail**.

The held-out run confirms the verdict: RT width selection stays off by default.

CUR (the current default with `library12_full.pqp`, self-calibrating) on all three raws:
- FDP_c ≤ 0.0096 on all 15 runs.
- Medians: 1804 / 1909 / 1907.
- These are 87 % / 83 % of DIA-NN's 2067 / 2301 on raws 01/02.

## Follow-up (16 November 2026)
The user judged a post-hoc RT filter the wrong place for this: the scorer should learn it. So
`-rt_width_selection` was removed (revert e560b3817a), not kept as an option.

The next direction is an RT-dependent extraction window; see `PRECHECK_RTPROFILE.md`.
