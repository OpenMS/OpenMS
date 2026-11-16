# RESULT_ENSEMBLE: seed-ensemble scorer (design runs, verdict)

Written 16 November 2026, after the raw 01/02 runs and before any raw 03 run.
- **Plan:** `PLAN_ENSEMBLE.md` (sha256 file present, docs 1414a163de).
- **Code:** OpenMS 09aca5207c. `CrossFitScorer::Options::ensemble` and OpenDIAWorkflow `-ensemble`. The unit test shows
  m = 1 is bit-identical to the single model, and m = 3 is thread-independent.
- **Container:** `opendia-ens-09aca5207c.sif` (sha256 5db38ced…), on codon with 8 threads and a node-local HOME.
- **Arms:** CUR = `-ensemble 1`, NEW = `-ensemble 3`, as separate runs. Self-calibrating, arm-E `library.pqp` plus
  `library12_full.pqp`, seeds 0-4.
- **Check:** the CUR runs reproduce the PLAN_RTWIN2 CUR runs exactly (same q ≤ 0.01 counts and N on all 10).
- **Evaluation:** `eval_run.py`. The tool's q equals the recompute on every file.

## N@FDP_c 1 % (FDP_c at q ≤ 0.01)
| run | arm | s0 | s1 | s2 | s3 | s4 | median | range |
|---|---|---|---|---|---|---|---|---|
| raw01 | CUR | 1815 (.0079) | 1804 (.0082) | 1784 (.0047) | 1728 (.0095) | 1827 (.0046) | **1804** | 99 |
| raw01 | NEW | 1817 (.0079) | 1807 (.0092) | 1822 (.0047) | 1808 (.0092) | 1805 (.0091) | **1808** | **17** |
| raw02 | CUR | 1894 (.0088) | 1928 (.0044) | 1935 (.0033) | 1867 (.0046) | 1909 (.0012) | **1909** | 68 |
| raw02 | NEW | 1891 (.0046) | 1920 (.0045) | 1933 (.0023) | 1878 (.0046) | 1922 (.0034) | **1920** | 55 |

Medians at the other FDP_c levels:

| run | arm | 0.5 % | 2 % |
|---|---|---|---|
| raw01 | CUR | 1616 | 1993 |
| raw01 | NEW | 1684 | 2012 |
| raw02 | CUR | 1800 | 2018 |
| raw02 | NEW | 1816 | 2035 |

## Criteria
| criterion | raw01 | raw02 |
|---|---|---|
| 1. NEW median ≥ CUR median + 1 % | 1808 vs 1822 needed (+0.2 %): **fail** | 1920 vs 1928 needed (+0.6 %): **fail** |
| 2. NEW FDP_c ≤ 0.010 every seed | max 0.0092: pass | max 0.0046: pass |
| 3. NEW seed range ≤ CUR range | 17 ≤ 99: pass | 55 ≤ 68: pass |

**Verdict: FAIL** on criterion 1. Per the plan, `-ensemble` is removed and nothing is retuned.

## Reading
- **The ensemble mostly removes bad seeds rather than raising the typical one.** On raw01 the seed range falls from 99
  to 17, and seed 3 goes 1728 → 1808. The median moves only +4 / +11.
- **The gain is larger away from the 1 % point:** +68 / +16 at 0.5 % and +19 / +17 at 2 %.
- **Offline probe against the real runs:** the probe predicted +26 / +34 at 1 %. Part of that came from comparing
  ensembles against a single-seed median that included the arrays' bad seeds.
- **Cost:** scoring time 17-22 s → 47-81 s; total run time is within run-to-run noise (343-509 s against 319-453 s).
  Memory is unchanged (21.9 GB).
- **For later:** the stability benefit is real. If run-to-run reproducibility becomes an acceptance requirement, this is
  the cheapest way to get it, but it would need its own plan.

## Raw 03 (held out, run once after the verdict was committed, docs 1b6388a27e)
CUR is the raw 03 CUR of PLAN_RTWIN2, which is the same code path as `-ensemble 1`.

| arm | s0 | s1 | s2 | s3 | s4 | median | range |
|---|---|---|---|---|---|---|---|
| CUR | 1907 (.0096) | 1831 (.0058) | 1864 (.0068) | 1913 (.0084) | 1957 (.0076) | **1907** | 126 |
| NEW | 1832 (**.0106**) | 1912 (.0088) | 1897 (.0087) | 1938 (.0096) | 2011 (.0097) | **1912** | 179 |

- Median +0.3 %, and seed 0 has FDP_c 0.0106. Raw 03 shows neither the gain nor the stability.
- The verdict stands: removed (revert bc30ae9394).
