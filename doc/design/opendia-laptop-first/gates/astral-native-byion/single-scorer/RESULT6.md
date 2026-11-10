# Round 6 result: transition refinement in the native CNN + LDA. Both arms FAIL

Evaluated with `evaluate6.py` (`EVALUATION6.txt`) under PLAN6.md (sha256 in `PLAN6.md.sha256`; pre-run column correction in `PLAN6.amended.sha256`). HeLa default equivalence: `crossfit.tsv` IDENTICAL (`hela-default-cmp6`). All q-derived counts are decoy-count estimates.

| Arm | Raw | N at FDP_c 0.5 / 1 / 1.5 / 2 % | q ≤ 0.01: N, FDP_c | q ≤ 0.05: N, FDP_c | Rule |
|---|---|---|---|---|---|
| base | 01 | 899 / 1029 / 1270 / 1321 | 1068, 0.013 | 1638, 0.065 | (reference) |
| base | 02 | 1230 / 1298 / 1331 / 1372 | 1236, 0.007 | 1788, 0.074 | (reference) |
| R1 (+refine) | 01 | 914 / **1028** / 1138 / 1261 | 1016, 0.008 | 1585, 0.076 | fail (rules 1, 3) |
| R1 (+refine) | 02 | 1055 / **1160** / 1291 / 1317 | 1148, 0.009 | 1724, 0.067 | fail (rule 1) |
| R2 (+mask 0.75) | 01 | 473 / **520** / 580 / 716 | 617, 0.016 | 1011, 0.072 | fail (rules 1–3) |
| R2 (+mask 0.75) | 02 | 482 / **679** / 727 / 757 | 632, 0.006 | 957, 0.057 | fail (rule 1) |

Each run took about 13 min with a peak of 23.8 GB.

## Internal check (R1): a planning error

PLAN6 copied the comparison over crossfit columns 1–7 and 9. Column 7 (`precursor_best`) is chosen by the stage-2 score, so it differs by design. The earlier byion CMP already noted this. The planned check therefore fails for both raws. A post-hoc check of columns 1–6 and 9 (sample, ids, decoy, fold, logit, window) and of feature columns 1–16 is IDENTICAL on both raws (`CMP6.txt`). This means the CNN and the old features were unchanged. The verdict does not depend on this: R1 fails rule 1 on both raws either way.

## Mechanism (reported only)

AUC of entrapment winners over decoy winners, over in-window winners by `stage2_score`:

| Column | Raw 01 R1 | Raw 02 R1 | R2 (both raws) |
|---|---|---|---|
| refine_quantitative | 0.508 | 0.507 | 0.503 |
| refine_best_reference_r | 0.512 | 0.512 | 0.503 |
| transition_support | 0.525 | 0.524 | 0.524–0.523 |
| cnn_z | 0.516 | 0.516 | 0.515 / 0.514 (base 0.518 / 0.519) |

- **Part 1 is supported.** On both raws the refine counts separate false targets from decoys about 3× less than `transition_support` (|AUC − 0.5| is 0.007–0.012 against 0.024–0.025).
- **Part 2 is supported only weakly.** With the mask, `cnn_z` moves 0.003–0.005 toward 0.5.

The interference reading holds partly. Refined evidence is less "target-like for absent peptides", but it is also weaker evidence for present peptides.

## Why the mask fails

At r < 0.75 the mask zeroed 8.6M channels in 2.53M of 2.82M windows, about 3 of 6 fragments per window. With only six fragments, and a peak of a few 0.94-s scan bins, the median consensus is noisy. The cut then removes real signal on most windows. That includes true peptides, which is why IDs halve. EncyclopeDIA applies the cut to quantification after detection, with many more transitions. It does not apply it to detection.

## Why R1 loses IDs

The six refine columns are strongly collinear, and the LDA gives them large opposite-signed weights. For example, in R2 fold 0, `refine_mean_r` is −122 and `refine_best_reference_r` is +113. R1 is more conservative at q ≤ 0.01 (FDP_c 0.008 / 0.009), and better at FDP_c 0.5 % on raw 01. It loses IDs at 1–2 % on raw 02.

## Not concluded

- A milder cut, such as down-weighting fragments instead of zeroing them, or a cut only when r < 0.5.
- A single refine column instead of six.
- Refinement over a wider window. None of these were tested.
