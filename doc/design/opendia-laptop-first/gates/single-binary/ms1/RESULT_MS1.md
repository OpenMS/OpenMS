# RESULT_MS1: MS1 isotope evidence block (design runs, verdict)

Written 16 November 2026, after the raw 01/02 runs and before any raw 03 run.
- **Plan:** `PLAN_MS1.md` (sha256 file present, docs 933869a919).
- **Code:** OpenMS 51fb7cc7ec.
  - `OpenDIAMS1`: envelope and 5 features.
  - `TraceCollector::enableIsotopeTraces`.
  - OpenDIAWorkflow `-ms1_features`.
  - `OpenDIAMS1_test` checks the envelope and the features against the Python precheck on 5 real arm-E windows;
    24/24 ctests pass.
- **Container:** `opendia-ms1-51fb7cc7ec.sif` (sha256 27293e32…), on codon with 8 threads and a node-local HOME.
- **Arms:** CUR = no flag (24 inputs), NEW = `-ms1_features` (29 inputs). Both are separate runs: self-calibrating,
  arm-E `library.pqp` + `library12_full.pqp`, seeds 0-4.
- **Checks:**
  - CUR reproduces the earlier CUR runs exactly.
  - The tool's q equals the recompute on every file.

## N@FDP_c 1 % (FDP_c at q ≤ 0.01) and IDs at q ≤ 0.01
| run | arm | s0 | s1 | s2 | s3 | s4 | median N@1 % | median q ≤ 0.01 |
|---|---|---|---|---|---|---|---|---|
| raw01 | CUR | 1815 (.0079) | 1804 (.0082) | 1784 (.0047) | 1728 (.0095) | 1827 (.0046) | **1804** | 1709 |
| raw01 | NEW | 1883 (.0049) | 1858 (.0047) | 1920 (.0057) | 1919 (.0046) | 1907 (.0047) | **1907 (+5.7 %)** | 1710 (+0.1 %) |
| raw02 | CUR | 1894 (.0088) | 1928 (.0044) | 1935 (.0033) | 1867 (.0046) | 1909 (.0012) | **1909** | 1811 |
| raw02 | NEW | 1919 (.0054) | 2012 (.0022) | 1961 (.0063) | 1920 (.0085) | 1949 (.0044) | **1949 (+2.1 %)** | 1861 (+2.8 %) |

N@FDP_c 0.5 % and 2 % medians:

| run | arm | 0.5 % | 2 % |
|---|---|---|---|
| raw01 | CUR | 1616 | 1993 |
| raw01 | NEW | 1727 | 2003 |
| raw02 | CUR | 1800 | 2018 |
| raw02 | NEW | 1846 | 2077 |

- Every NEW seed is above every CUR seed on raw01.
- The seed range is 99 → 62 on raw01 and 68 → 93 on raw02.

## Criteria
| criterion | raw01 | raw02 |
|---|---|---|
| 1. NEW median N@1 % ≥ CUR − 1 % | 1907 vs 1786: pass | 1949 vs 1890: pass |
| 2. NEW FDP_c ≤ 0.010 every seed | max 0.0057: pass | max 0.0085: pass |
| 3. NEW median q ≤ 0.01 ≥ CUR + 2 % | 1710 vs 1743 needed (+0.1 %): **fail** | 1861 vs 1847: pass |
| 4. AUC(E > ED) NEW ≤ CUR + 0.02 (`CRIT4.txt`) | 0.5307 vs 0.5274 allowed (CUR 0.5074): **fail** | 0.5086 vs 0.5393: pass |

**Verdict: FAIL** on criteria 3 and 4 on raw01. Per the plan, `-ms1_features` is removed and nothing is retuned.

## Reading (does not change the verdict)
- **The plan expected the wrong effect.** The offline probe showed the gain in the q ≤ 0.01 count; the full workflow
  shows it in the ranking at equal true FDR instead: +103 / +40 at FDP_c 1 %, the largest gain since block I.
  - On raw01 the q-values stay conservative: FDP_c 0.005 at q ≤ 0.01 against 0.008 for CUR.
  - So the count at q ≤ 0.01 does not rise, which fails criterion 3.
- **Criterion 4 fails narrowly on raw01:** 0.531 against 0.527.
  - Entrapments in the rejected band rank slightly higher under NEW. This is the benchmark caveat in the plan:
    entrapments share their parent's m/z and envelope.
  - Over all winners below the cut, the AUC(entrapment > decoy) is unchanged (0.500 / 0.500).
  - FDP_c at q ≤ 0.01 is lower under NEW on raw01.
  - The leak is small, but it is a leak.
- **The C++ block differs from the precheck in one input.** The C++ fragment-sum trace comes from the 6-fragment tensor
  (0.02 Da windows), as the plan specifies. The precheck used `xic.npy` (±17.5 ppm around the calibrated offset).
  - Across 919k matched raw01 windows, m1_corr_xic agrees at Pearson 0.73 and m1_apex_dt at 0.94.
  - The MS1-only features agree at 0.9997 or better.
  - The extraction, the isotope traces and the envelope are therefore exact. Only the fragment reference differs.
- **Cost:** 5 more features. Extraction time is within run-to-run noise (190-430 s against 148-324 s), and memory is
  unchanged (21.9 GB).
- **Option for the user, not taken here.** The goal metric (equal-FDR N) improved on both raws. A new plan with
  criteria written for that effect would need fresh data for its verdict: raw 03 below, or other runs. That is a
  decision for the user, not a retuning of this plan.

Raw 03 (held out) follows below, reported descriptively.

## Raw 03 (held out, run once after the verdict was committed, docs 05600821f6)
CUR is the raw 03 CUR of the earlier plans, which uses the same 24-input code path.

| arm | s0 | s1 | s2 | s3 | s4 | median N@1 % | median q ≤ 0.01 |
|---|---|---|---|---|---|---|---|
| CUR | 1907 (.0096) | 1831 (.0058) | 1864 (.0068) | 1913 (.0084) | 1957 (.0076) | **1907** | 1849 |
| NEW | 1898 (.0055) | 1835 (**.0100**) | 1951 (.0077) | 1941 (.0066) | 1893 (.0077) | **1898 (−0.5 %)** | 1811 (−2.1 %) |

- **Raw 03 does not confirm the gain.** N@1 % changes by −0.5 % and the q ≤ 0.01 count by −2.1 %. One seed sits at
  FDP_c 0.0100, exactly at the limit.
- The raw01 / raw02 gain (+5.7 % / +2.1 %) therefore does not carry over to the held-out run.
- The verdict stands: removed (revert below). The block's code and tests stay in git history (51fb7cc7ec).
