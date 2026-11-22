# RESULT_GAP2: the gap after sibling rescoring (22 Nov 2026)

The gate is `PRECHECK_GAP2.md`.
- **Gate timing:** committed locally (`a9119bb7ab`, 05:47:49Z) before the first job started (05:50:55Z). The push failed three times on network timeouts and landed at about 05:56Z, while the jobs were running.
- **Runs:** all 12 runs and 4 analyses completed (`gap2-20261122`). Full output is in `EVAL_GAP2.txt` (part A) and `EVAL_GAP2_B.txt` (part B).
- **Features:** the default tables now carry 49 features.

## Part A: arms on the new default (seed 0)

HeLa IDs at own q ≤ 0.01. At r7 each cell is followed by N @ FDP_c 1 % in brackets.

| Cell | DIA-NN | A tool | B label-free GBM | Bt hard decoys | C DIA-NN ceiling |
|---|---:|---:|---:|---:|---:|
| D r0 | 27,901 | 25,863 (0.93) | 28,078 (+8.6 %) | 27,919 (+7.9 %) | 24,351 (−5.8 %) |
| H r0 | 44,866 | 39,365 (0.88) | 41,904 (+6.4 %) | 41,470 (+5.3 %) | 40,850 (+3.8 %) |
| D r7 | 25,881 | 19,932 [19,269] | 20,718 [20,030] | 22,175 [**21,694**] | 16,473 [16,205] |
| H r7 | 41,231 | 35,304 [34,953] | 36,799 [35,936] | 38,091 [**37,534**] | 34,053 [33,363] |

**Bounds at r7 (lower / combined):**

| Run | A | B | Bt |
|---|---|---|---|
| D | 1.17 / 1.34 % | 1.20 / 1.37 % | 1.15 / 1.31 % |
| H | 1.16 / 1.33 % | 1.23 / 1.40 % | 1.29 / 1.47 % |

**Funnel of DIA-NN-only HeLa IDs:**

| Cell | Missed | A_low | B2 | C |
|---|---:|---:|---:|---:|
| D r0 | 10,175 | 61.6 % | 10.8 % | 27.6 % |
| H r0 | 11,705 | 37.6 % | 12.2 % | 50.3 % |
| D r7 | 12,096 | 64.0 % | 11.4 % | 24.6 % |
| H r7 | 12,018 | 49.1 % | 13.0 % | 37.9 % |

## Part B: candidate arms against the new default (seed 0)

| Arm | D r0 HeLa | H r0 HeLa | D r7 N @ FDP_c 1 % | H r7 N @ FDP_c 1 % | Lower bound change (r7) | Wall (max) | Rule B |
|---|---:|---:|---:|---:|---|---|---|
| V7M6 | +16.0 % | +3.3 % | +13.3 % | +3.6 % | −0.07 / −0.03 pp | 1.52× (D r7) | fail (wall) |
| V8M8 | +23.7 % | +8.8 % | +2.9 % | +6.3 % | −0.12 / −0.06 pp | 2.17× (D r7) | fail (wall) |

## Reading

1. **Direction rule:** branch 1 applies.
   - G_Bt at r7 N @ FDP_c 1 % is **D +12.6 %** (19,269 to 21,694) and **H +7.4 %** (34,953 to 37,534). Both are ≥ +5 %, so the next C++ work is a scorer that can use hard negatives.
   - On H, the Bt lower bound rises by +0.13 pp, which is within any +0.2 pp bar but is recorded.
2. **Siblings moved the gap.**
   - **D r0:** the tool now reaches 0.93 of DIA-NN, up from 0.68.
   - **DIA-NN ceiling is below the tool:** this holds at D r7, H r7 and D r0. Training on DIA-NN labels no longer helps, because the tool's own sibling-aware labels are at least as good.
   - **What remains:**
     - **A_low (D):** right peaks scored too low, 64 % of D's missed IDs at r7.
     - **No candidate, C (H):** 38-50 % of H's missed IDs.
3. **Candidate arms.** V7M6 now gains much more on D (r7 +13.3 %, against +3.5 % on the old default). It misses the wall bar on one cell only (1.52× against 1.5×; it was 0.79× on H r7, so node noise is likely). By the rules it fails, and I do not move the bar. A cheaper split, or re-measuring wall time on matched nodes, needs your decision.
