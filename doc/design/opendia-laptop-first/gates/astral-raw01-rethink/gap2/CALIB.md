# Why decoy-estimated q is anti-conservative (raws 01 and 02, OSW-like scorer)

This is a design diagnostic. Entrapment is used here only to diagnose. Any choice made from it must be confirmed on an untouched run (raw 03).

The scripts are `calib.py`, `calib_curve.py`, `calib_tail.py`, `calib_borrow.py`, `calib_hits.py`, `calib_il.py`, `calib_shift.py`, `calib_struct.py`, `calib_entr.py`, `calib_entr2.py`, `calib_lowmass.py`, `calib_y1.py`, `calib_scan.py`, `calib_mitig.py` and `matched_fdp02.py`. Each writes a matching `.txt` output.

## 1. Size and shape (`calib_curve.txt`)

Decoy-estimated FDR against entrapment FDP_c along the accepted list. The 90% intervals are exact Poisson.

| run, scorer | q ≤ 0.01: N, FDP_c [90% CI] | q ≤ 0.05: N, FDP_c [90% CI] |
|---|---|---|
| 01 LDA | 993, 0.010 [0.004, 0.021] | 1429, 0.058 [0.044, 0.075] |
| 01 MLP | 1222, 0.012 [0.005, 0.022] | 1691, 0.078 [0.062, 0.095] |
| 02 LDA | 1115, 0.024 [0.014, 0.037] | 1463, 0.079 [0.062, 0.098] |
| 02 MLP | 1231, 0.020 [0.011, 0.032] | 1573, 0.072 [0.057, 0.090] |

**The excess is systematic, not noise.** At q ≤ 0.05 every interval lies above 0.05 on both runs, for both learners. At q ≤ 0.01 it is significant on raw 02 only. OSW+Percolator's own m_score ≤ 0.01 export has FDP_c 0.032 (raw 01) and 0.034 (raw 02) on the same library and decoys (`matched_fdp02.txt`).

**Bound caveat.** FDP_c is the *combined* entrapment estimate, which is an estimated upper bound (Wen et al., Nat Methods 2025). With r ≈ 1 it is about twice the *lower bound* N_e/N. The paired bound cannot be used here: entrapment peptides are not shuffles of this library's human peptides, and only 4.2% share a human precursor's composition. Lower bounds at q ≤ 0.01:

| run | LDA | MLP | OSW export |
|---|---|---|---|
| raw 01 | 5/993 = 0.005 | 7/1222 = 0.006 | 23/1462 = 0.016 |
| raw 02 | 13/1115 = 0.012 | 12/1231 = 0.010 | 25/1491 = 0.017 |

At 1% our scorer is therefore *inconclusive* in Wen et al.'s sense: the upper bound is above the line and the lower bound is at or below it. OSW's export is invalid even by the lower bound. DIA-NN 2.3.1, with its own different library, reports FDP_c 0.007 at 1% on all three runs (lower bound 0.0036). So valid control is achievable on these raw files. The direct evidence of a real excess comes from section 2, not from FDP_c > 1% alone.

## 2. Where it comes from (`calib_tail.txt`)

Before competition, let S_e be the share of entrapment target winners above a score cut, and S_d the same share for decoy winners. Valid target-decoy competition needs S_e/S_d ≤ about 1. The 90% Poisson band is in brackets.

| run, scorer | cut at q ≤ 0.01 | cut at q ≤ 0.02 | cut at q ≤ 0.05 |
|---|---|---|---|
| 01 LDA | 1.22 [0.32, 2.21] | 1.12 | 1.12 [0.83, 1.40] |
| 01 MLP | 1.14 | 1.67 [1.10, 2.25] | 1.52 [1.21, 1.83] |
| 02 LDA | 2.55 [1.38, 3.75] | 1.60 | 1.55 [1.21, 1.89] |
| 02 MLP | 2.14 [1.12, 3.19] | 1.65 | 1.41 [1.10, 1.72] |

**In the score tail, false targets beat decoys by about 1.1–1.6×, and by up to 2.5× near the 1% cut on raw 02.** The LDA shows it too, so this is mainly not MLP overfitting. It is a mismatch between the decoy model and false targets. The MLP adds to it at the 2–5% cuts on raw 01.

## 3. What it is not

- **Not a single feature.** Over all winners, every subscore has an entrapment-vs-decoy AUC within 0.48–0.53 (`calib_auc_raw0*.tsv`). Among the winners above the 5% cut, no subscore differs at p < 0.02 (`calib_hits.txt`, n = 41–85). The closest is MS1 shape on raw 02 (p 0.02–0.08), and MS1 adds nothing to N. So dropping a feature is not a supported fix.
- **Not true peptides in disguise.** No library entrapment sequence equals a human one under I = L or Q/K, N/D swaps (`calib_il.txt`). Only 1 of the tail hits on raw 02 was reported by DIA-NN at 1%.
- **Borrowed signal: not shown on raw 01, weakly suggested on raw 02.** Here neighbour overlap means the number of top-6 fragments shared with a human precursor within ±1 Th and ±30 s (`calib_borrow.txt`). On raw 02 the MLP's 1% entrapment hits share ≥ 2 fragments with a co-eluting human precursor in 2 of 12 cases, and the 5% hits in 6 of 56. The decoys above the same cuts show 0 of 11 and 2 of 78. Across random precursors the rates are 6.0% (entrapment) and 4.9% (decoys). On raw 01 there is no enrichment: 0 of 5 and 0 of 41 for the LDA. With these counts this is a hypothesis, not a finding.

## 4. Mitigation test (predeclared, `calib_mitig.txt`)

**Pass rule:** q ≤ 0.01 FDP_c ≤ 0.015 on **both** runs, with N at FDP_c 1% within 3% of the single-seed MLP or better.

| variant | raw 01: q ≤ 0.01 N (FDP_c), N at FDP_c 1% | raw 02: q ≤ 0.01 N (FDP_c), N at FDP_c 1% |
|---|---|---|
| A LDA | 993 (0.010), 1011 | 1115 (0.023), 885 |
| B MLP seed 0 | 1222 (0.012), 1221 | 1231 (0.020), 1040 |
| C MLP 5-seed ensemble | 1171 (0.010), **1234** | 1251 (0.024), **1136** |
| D rank average of A and C | 1138 (0.012), 1099 | 1208 (0.027), 1009 |

**The test fails.** The ensemble is the better ranker: N at matched true error goes up on both runs (1234 against 1221, and 1136 against 1040), and it calibrates raw 01 better at 2–5%. It does not fix raw 02 at 1% (FDP_c 0.024). On raw 02 the ensemble still reaches only 88% of OSW+Percolator's 1287.

## 5. Structural checks (`calib_struct.txt`, `calib_entr*.txt`, `calib_lowmass.txt`, `calib_y1.txt`, `calib_scan.txt`)

- **Groups match on precursor m/z, charge, length, library RT, oxidation, C content, fragment type and ordinal, XIC points and RT span.**
- **An out-of-scan-range artefact in the library.** No signal exists below 150.06 Th in these XICs, which is the acquisition's lower limit. Yet the library's top-6 fragments include ions below 150 Th:
  - in 9.7% of decoys, against 5.7% of entrapment and 5.2% of human precursors;
  - the decoy set even includes transitions at 55–100 Th, below the 100 Th assay limit;
  - the y1 of a K-terminal peptide is 147.11 Th, and pseudo-reversal moves the decoy's terminal residue around, so more decoys carry a dead y1.

  A decoy with a dead slot never reaches the score tail (0 of 23,446 above the 5% cut), which explains why decoys show more empty traces. **It is a real defect, but it does not explain the excess.** Dropping every pair with a slot below 150 Th changes S_e/S_d at the 1%/5% cuts only slightly: LDA raw 01 1.22/1.12 → 1.22/1.07, LDA raw 02 2.55/1.55 → 2.82/1.56, ensemble raw 02 2.67/1.22 → 2.66/1.20. FDP_c at q ≤ 0.01 is unchanged.
  - **Fix for the native pipeline:** pick the top-N fragments only within the acquired scan range, and apply the same rule to targets and decoys. This also returns a dead slot to 5.2% of human targets, a sensitivity gain.
- **Entrapment traces are less often empty than both human and decoy traces at the same fragment m/z** (for example, 300–500 Th: entrapment 0.069, human 0.099, decoy 0.091). Re-weighting entrapment winners to the decoys' strata (empty-trace count × candidate count) lowers S_e/S_d by only about 5% (`calib_entr2.txt`). Entrapment is therefore slightly more signal-prone than human false targets, which biases FDP_c a little upward, but not enough to explain the excess.

## Reading

1. The anti-conservativeness is a property of **decoy winners under-representing false targets in the score tail**, by about 1.1–1.6×. The pipeline's ranker only adds to it. OSW+Percolator, with the same decoys, is worse: invalid even by the lower bound. No ranker change alone has fixed it. At 1% our scorer sits in the inconclusive zone between the entrapment bounds (section 1). The cause is still unidentified. These have been checked and are not it: single features, I/L, the precursor mass shift, borrowed signal (weak), dead low-m/z slots, and entrapment stratum mismatch (each about 5% at most).
2. **The decoy precursor mass shift is not supported as the cause.** The library's decoys are pseudo-reversed with a precursor mass shift of ±12, ±28 or ±44 Da (99.9% of pairs; checked in the PQP). The entrapment peptides carry no shift; they are scored at their own m/z. If the shift caused the excess, tail survival should vary with the sign or size of the shift (`calib_shift.py`). It does not vary consistently:
   - for negative against positive shifts, S_e/S_d at the 5% cut is 1.13 against 1.10 (raw 01 LDA) and 1.35 against 1.83 (raw 02 LDA);
   - the |12| and |44| groups hold only 0–7 decoys above the cut.

   The hypothesis is not excluded, but there is no in-data support. A clean test would need regenerated decoys at the target's own m/z.
3. **Levers, in order of how directly they hit the cause:**
   - (a) A **decoy-generation test**: regenerate decoys without the mass shift (and, separately, as shuffles instead of pseudo-reverse), extract and rescore, then repeat the tail test. This needs a new library and an extraction run on raws 01 and 02. Native OpenDIA extraction is the natural vehicle, since OSW's XICs are fixed to these decoys. The in-data shift test above found no support, so the expected effect is uncertain.
   - (b) **More decoys per target** (for example two decoy families, with q computed from their mean). This lowers the variance of the tail estimate but does not remove a bias.
   - (c) Report **a decoy-count estimate corrected by a measured tail factor** only as a diagnostic, never in production. Entrapment-based calibration would not be available on real data.
   - (d) Keep the **5-seed ensemble** as the ranker. It is the best N at matched true error on both runs.
4. The pass rules stay as they are: q ≤ 0.01 must hold at FDP_c ≤ 0.015, and on raw 02 nothing so far meets that.
