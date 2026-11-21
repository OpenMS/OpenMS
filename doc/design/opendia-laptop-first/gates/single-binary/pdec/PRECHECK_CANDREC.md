# PRECHECK_CANDREC: why do DIA-NN's peaks get no OpenDIA candidate, and which fix recovers them? (21 Nov 2026, diagnostic)

**Question.** RESULT_CEIL0 found that 27 % (D) and 48 % (H) of DIA-NN-only HeLa IDs at r0 have no OpenDIA window within 6 s of DIA-NN's apex. That share is flat across intensity quartiles. PRECHECK_NOCAND (older code, r7 D) traced 88 % of such misses to the search valley rule merging the right maximum into a neighbour. Widening the picker recovered the candidates, but only +1.8 % IDs at 2× wall time. Since then the RT window, the RT refinement and H have all changed. This precheck re-measures the causes on the current default and asks whether a candidate fix is worth a plan.

## Part A: causes (no new runs)

- **Data:** the four PRECHECK_CEIL0 `-features_out` runs (D/H × r0/r7, seed 0, byte-identical to the default) and the DIA-NN HeLa q ≤ 0.01 IDs from the same library level.
- **Script:** `precheck_candrec_a.py`. The refined RT prediction is recovered per precursor from `rt_apex ± abs_library_rt_error`. The refined tolerance is read from the run log.
- **Classes,** for C IDs (missed by the tool and no window within 6 s):
  - **rt_out:** DIA-NN's RT is outside the refined prediction ± tolerance.
  - **near:** the nearest window is 6-15 s away.
  - **cap:** inside the window, 5 windows already, none near.
  - **picker:** inside the window, fewer than 5 windows, none near.

## Part B: candidate fixes (end to end, only after part A)

- **Arms on the r0 and r7 libraries, D and H, seed 0:** the default against the existing advanced options, chosen by the dominant part-A class:
  - picker or cap dominant: `-valley_fraction 0.8 -max_peaks 8`.
  - rt_out dominant: a wider RT tolerance, which needs a new option.
- **Decision rule (written before part A):**
  1. **Promote to a frozen plan** if the best arm gives ≥ +3 % HeLa IDs at own q ≤ 0.01 at r0 on both runs, and ≥ +2 % N@FDP_c 1 % at r7 on both runs. The r7 lower bound may rise by no more than +0.2 pp, and wall time must stay ≤ 1.5× the default's.
  2. **Record and stop** if the arm recovers right candidates for ≥ 70 % of part-A C IDs but fails rule 1. The bottleneck is then evidence, not candidates, and work moves to A_low features and the Bt scorer port.
  3. **The fix does not work** if it recovers right candidates for < 50 % of C IDs. Part A's attribution is then wrong and must be revisited before any further candidate work.

No code or default changes follow from part A.
