# RESULT_JOINT2: hard decoys plus sibling evidence gives +20-26 % at true FDP 1 % on every run (21 Nov 2026, diagnostic)

The gate is in RESULT_JOINT, section "Follow-up fixed before running" (`62d6becd5f`). The script is `precheck_joint2.py` (sha256 `d410c287…`). The data is the RESULT_BT156 r ≈ 1.56 feature tables (D/H seeds 0-1, r = 1.1225) and the PRECHECK_CEIL0 r7 tables (seed 0). Labels are used only for evaluation. The full output is in `EVAL_JOINT2.txt`.

N at combined FDP 1 % (lower bound at own q ≤ 0.01). Gains are against A, the tool score.

| Run | A tool | B | Bs (+ siblings) | Bt (hard decoys) | **Bts (hard decoys + siblings)** | Bts gain | Lower Δ |
|---|---:|---:|---:|---:|---:|---:|---:|
| D r≈1.56 s0 | 16,678 (0.89 %) | 17,265 | 17,223 | 17,434 | **20,152 (0.90 %)** | **+20.8 %** | +0.01 pp |
| D r≈1.56 s1 | 17,038 (0.88 %) | 17,271 | 17,638 | 17,524 | **21,494 (0.88 %)** | **+26.2 %** | 0.00 pp |
| H r≈1.56 s0 | 32,310 (0.66 %) | 33,711 | 38,248 | 34,378 | **38,958 (0.76 %)** | **+20.6 %** | +0.10 pp |
| H r≈1.56 s1 | 32,415 (0.68 %) | 33,674 | 38,406 | 34,244 | **38,900 (0.76 %)** | **+20.0 %** | +0.08 pp |
| D r7 s0 | 15,465 (1.12 %) | 15,387 | 14,780 | 16,131 | **18,672 (1.27 %)** | **+20.7 %** | +0.15 pp |
| H r7 s0 | 30,543 (1.20 %) | 31,250 | 36,588 | 32,353 | **37,500 (1.24 %)** | **+22.8 %** | +0.04 pp |

**Rule:**

| Criterion | Result | |
|---|---|---|
| ≥ +3 % on all 4 r≈1.56 runs | +20.0 to +26.2 % | yes |
| Mean lower-bound rise ≤ +0.2 pp per run | D +0.005, H +0.09 pp | yes |
| ≥ +2 % at r7 on D and H | +20.7 %, +22.8 % | yes |

**Verdict: PASS for Bts.** It goes to a C++ opt-in.

## Reading

- **The two changes need each other on D.**
  - Siblings alone (Bs) add nothing on D: −0.3 / +3.5 % at r≈1.56 and −4.4 % at r7.
  - Hard decoys alone (Bt) add +2.9 to +4.5 %.
  - Together they add +21-26 %. With all decoys as negatives, the learner spends its capacity separating the bulk of decoys and never uses sibling support at the cut. Against the hard decoys only, sibling co-elution is exactly what separates a faint true peak from a target-like decoy.
  - On H the siblings carry most of the gain alone (Bs +18 %), and hard decoys add about 2 points more.
- **Where this puts us against DIA-NN at r7, true FDP 1 %:**

  | Run | Bts | DIA-NN | Now |
  |---|---|---|---|
  | D | 18,672 = **72 %** of 26,006 | | 59 % |
  | H | 37,500 = **91 %** of 41,272 | | 74 % |

  H would pass the 90 % target. D would close about a third of its remaining gap.
- **Calibration is checked by entrapment, not assumed.**
  - Sibling evidence is symmetric. Decoy peptides carry decoy precursors at the same charges as their targets. Entrapment peptides carry entrapment siblings, so chance co-elution boosts entrapments exactly as it would boost false targets.
  - N@FDP_c 1 % is the entrapment-controlled count. The decoy-based lower bound rises at most +0.15 pp (D r7), within the +0.2 pp limit.
- **Why siblings carry new information.** A true peptide seen at two charge states elutes at the same time in both. One faint precursor is then backed by its bright sibling, which is evidence that no per-precursor feature can contain. The DIA-NN-labelled ceiling with siblings rose +21.5 % on D r0 (RESULT_JOINT), against less than +1 % for every per-precursor feature tested before.
- **Caveat on cross-fitting.** `s_best` uses the tool score of the sibling's window. A sibling may sit in a different library fold, so its score can come from a model trained on folds that include this precursor. That passes only target/decoy information, which every arm already uses, never entrapment or DIA-NN labels. It could still flatter the result slightly. The C++ port must compute sibling features from scores that are out-of-fold for the precursor, or assign siblings to one fold. Item (b) below measures the effect.

## What the C++ port has to get right

1. **Sibling features need a first-pass score.** The workflow already runs two passes (RT refinement). The sibling features (`s_n`, `s_co`, `s_best`, `s_best_d`) can be computed from the first pass's per-window scores and added to the second pass's scorer input. RESULT_JOINT2 used the final tool score. A first-pass score is weaker, so the gain may be smaller.
2. **The learner.** RESULT_BTCPP showed that the C++ network with hard decoys collapses its decoy whitening. Bts needs hard decoys. So either the sibling features have to work in the network with its default negatives (Bs, which failed on D with the GBM), or the scorer needs the redesign noted in RESULT_BTCPP: hard-decoy loss weighting with whitening on all decoys, or a GBM-like final stage.
3. **Recommended order:**
   - (a) A C++ opt-in for sibling features only, tested with the network at default negatives on r ≈ 1.56, two seeds. This is cheap and isolates point 1.
   - (b) A Python check of Bts with first-pass scores on the existing tables, to bound the loss from point 1.
   - (c) The scorer redesign for hard negatives.
