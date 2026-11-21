# PRECHECK_CEIL0: is the plain-library gap a scoring problem or a candidate and feature problem? (21 Nov 2026, diagnostic)

**Question.** RESULT_LIBSIZE showed that most of the gap is present on the plain library (r0): OpenDIA finds 0.683 (D) and 0.774 (H) as many HeLa IDs as DIA-NN. PRECHECK_SCORER36 found, on older code and at r7 only, that a learner trained on DIA-NN labels does no better than the tool. This precheck repeats that ceiling on the current default, at r0 and r7, on D and H. It answers two questions. Does our feature table carry the information needed for DIA-NN's extra IDs? Where in the funnel are those IDs lost?

## Setup (fixed before any result)

- **Runs:** current default (`opendia-pd10-rtrdef.sif`, sha256 `f334b5fc…`), seed 0, `-features_out`, on D and H with the r0 library (`libsize-20261121/r0/lib.pqp`) and the r7 library.
  - **Validity:** each `precursors.tsv` must be byte-identical to the matching seed-0 run without `-features_out`: `libsize-20261121/r0/opendia/{D,H}_s0` for r0, and `entlib-rtr-20261121/runs/{D,H}_rtr_s0` for r7. If it is not, the run is reported but not used.
- **Labels:** DIA-NN 2.2 HeLa IDs at q ≤ 0.01 from the same library level (`libsize-20261121/r0/diann/{D,H}` and `entlib-diann5-20261119/{D,H}`). They are used only for training the ceiling and for evaluation.
- **Script:** `precheck_ceil0.py` (sha256 in the result). The method is the same as `precheck_scorer36.py`.
- **Arms:**
  - **A:** the tool score.
  - **B:** a label-free GBM with 3 Percolator rounds.
  - **Bt:** B with hard decoys as negatives.
  - **C:** a GBM trained on DIA-NN labels. Positives are target windows within 6 s of DIA-NN's RT.
  - **Ct:** C with hard decoys as negatives.
- **Measure:** HeLa IDs at own q ≤ 0.01 under paired target-decoy competition. N at combined FDP 1 % is also reported at r7.
- **Funnel:** DIA-NN-only IDs are split into three classes:
  - **A_low:** the right window is the best window but scores too low.
  - **B2:** the right window exists but loses to another window.
  - **C:** there is no window within 6 s.
- **TPR:** of the right window at decoy FPR 1e-4, 1e-3 and 1e-2, by DIA-NN quantity quartile.

## Decision rule (written before any result)

The ceiling gain is G = max(C, Ct) / A − 1, in HeLa IDs at own q ≤ 0.01, at r0.

1. **The scorer is a real lever** if G ≥ +15 % on D or G ≥ +10 % on H. The features then hold information that our label-free scorer does not extract. The next step is a precheck of label-free scorer changes, starting with Bt and with features that DIA-NN labels mark as useful.
2. **The scorer is not the limit** if G < +5 % on both runs. Work then moves to the funnel. A_low drives new faint-peak evidence, and B2 and C drive candidate recovery (step 3). The largest class picks the next precheck.
3. Results between these bands are reported, and the funnel decides.
4. **Independent of rules 1-3:** if label-free B or Bt gives ≥ +5 % at r0 on both runs and ≥ +3 % at N@FDP_c 1 % at r7 on both runs, it becomes a candidate for its own r ≈ 1.56 precheck and frozen plan. Calibration at r7 must not worsen by more than +0.2 pp on the lower bound.

No code change follows from this precheck alone.
