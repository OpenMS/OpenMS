# PRECHECK_GAP2: where is the remaining gap after sibling rescoring? (22 Nov 2026)

**Why.** Charge-state sibling rescoring is now the default (OpenMS `e2e65938ce`). At r ≈ 7 it reaches 76 % of DIA-NN on D and 84 % on H at N @ FDP_c 1 %. The RESULT_CEIL0 funnel (A_low / B2 / C) and the candidate arms (RESULT_CANDREC, RESULT_SPLIT) were measured on the old default, so their shares and gains may have moved. This precheck re-measures both before choosing the next C++ work.

**Inputs.**
- **Container:** `opendia-pd10-sibdef.sif` (sha256 `78cc3a62…`, the new default). Seed 0 with `-features_out`.
- **Runs and levels:** D and H, at r0 (`libsize-20261121/r0/lib.pqp`) and r7 (`entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`).
- **DIA-NN references:** as in PRECHECK_CEIL0.

## Part A: diagnostic (no pass or fail)

- **Script:** `precheck_ceil0.py`, unchanged (sha256 `19955add…`), on the 4 default `-features_out` tables. The tables now carry 49 features, including the 4 sibling features.
- **Reported:**
  - **Arms:** A (tool), B (label-free GBM), Bt (label-free GBM with hard decoys), C (DIA-NN-labelled ceiling) and Ct.
  - **Funnel:** A_low / B2 / C of the DIA-NN-only HeLa IDs.
  - **TPR:** by quantity quartile.

## Part B: candidate arms re-baselined

- **Arms:**
  - **V7M6:** `-valley_fraction 0.7 -max_peaks 6`.
  - **V8M8:** `-valley_fraction 0.8 -max_peaks 8`.
- **Runs:** seed 0, same 4 cells, compared with the part A default runs. The evaluator is `precheck_gap2_b.py`, committed with this file. It reports HeLa IDs at own q ≤ 0.01, r7 N @ FDP_c 1 % with bounds, and wall time and peak RSS.
- **Rule B** (the PRECHECK_CANDREC rule 1 bars, unchanged). An arm passes if all hold:
  - r0 HeLa ≥ +3 % on D and H;
  - r7 N @ FDP_c 1 % ≥ +2 % on D and H;
  - r7 lower bound rises by ≤ +0.2 pp;
  - wall time ≤ 1.5× on every cell.
- **If an arm passes:** it goes to a frozen 5-seed r7 plan. If both pass, the one with the higher mean r7 gain goes.

## Direction rule (written before any result)

The next C++ work is chosen as follows. These are the measures used, from part A at r7:

- **G_Bt:** the Bt gain over A at N @ FDP_c 1 %.
- **G_C:** the C gain over A at HeLa own q ≤ 0.01.

The options, in order:

1. **Hard-negative scorer redesign.** Chosen if G_Bt ≥ +5 % on both D and H, so label-free scoring headroom remains.
2. **Candidate work.** Chosen if (1) does not hold and a part B arm passes rule B.
3. **Peak-level evidence.** Chosen if neither holds. This means a deconvolution precheck (A_low interference) or new features, and goes first if A_low is ≥ 50 % of the DIA-NN-only IDs on D or H at r7.

All numbers are reported whichever branch is taken.
