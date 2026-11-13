# RESULT_QUIRK: do decoys carry a tell that entrapments lack? (evaluation only, PAIRLIB run, D11FDWT)

Scripts: `scripts/decoy_quirk.py` (output `quirk/QUIRK_D11FDWT.txt`) and an inline rates script (output
`quirk/RATES_D11FDWT.txt`). Labels were used for evaluation only. No score, model or rule was changed.

## 1. At matched score, decoy and entrapment winners look the same (hypothesis as stated: not supported)
Population: winners above the q ≤ 0.20 cut. raw01: 727 decoy and 490 entrapment winners. raw02: 700 and 422.

| test | raw01 | raw02 |
|---|---|---|
| classifier on the 20 scorer inputs, decoy vs entrapment, CV AUC | 0.516 | 0.513 |
| same, score-matched | 0.505 | 0.504 |
| reference: human vs entrapment | 0.870 | 0.880 |
| largest single-feature score-matched AUC (boundary_library_correlation, f12_libcos) | 0.44–0.45 | 0.45–0.46 |

- **Above the cut the scorer's inputs cannot tell a decoy winner from an entrapment winner.** Two features lean the
  other way, mildly: entrapment winners fit the library pattern slightly *better* than decoy winners.
- So the scorer is not ranking entrapments above decoys at the same evidence.

## 2. Fewer decoys reach the tail in the first place: decoys match noise less often
Tail-entry rate per class = winners above a score cut / precursors of that class in the population. Classes:
- **H:** human target;
- **E:** entrapment;
- **D(H):** decoy of a human target;
- **D(E):** decoy of an entrapment.

Ratio of the entrapment rate to the mean decoy rate:

| cut | raw01 | raw02 |
|---|---|---|
| s ≥ 5.8 (≈ q 0.20) | 1.32 | 1.21 |
| s ≥ 7.2 (≈ q 0.10) | 1.41 | 1.32 |
| s ≥ 9.0 (≈ q 0.05) | 1.82 | 1.24 |

Null windows (|RT error| > 3σ, where no peptide is expected) show the same asymmetry, without any real signal
involved:

| class | mean support | support ≥ 4 | library cosine > 0.8 |
|---|---|---|---|
| H | 2.01 / 1.97 | 5.9 % / 5.4 % | 2.3 % / 2.3 % |
| E | 2.01 / 1.98 | 5.6 % / 5.4 % | 2.1 % / 2.1 % |
| D(H) | 1.93 / 1.89 | 4.9 % / 4.4 % | 1.7 % / 1.8 % |
| D(E) | 1.94 / 1.90 | 5.0 % / 4.6 % | 1.8 % / 1.6 % |

(raw01 / raw02 in each cell)

- On the scorer's inputs, a classifier separates decoy from target null windows with AUC 0.55 on both raws.
- **The decoy fragment sets hit random signal about 15–25 % less often than any target fragment set.** That includes
  entrapments, which are shuffles of the target sequences: shuffles behave exactly like real peptides.
- A high-scoring false target is mostly a lucky random match. Decoys are less lucky by construction, so they
  under-count false targets in the tail. This matches the measured rho of about 1.2–1.4 at q ≤ 0.10–0.20 directly.
  The tail of rho at q ≤ 0.05 (1.8 on raw01, 1.2 on raw02) is noisy.

## Reading
- The tell is in the **decoy library itself**, not in the scorer:
  - The decoys are pseudo-reversed with a K/R swap, are mass-shifted by about ±28 Da, and keep their target's
    intensities.
  - Their fragment m/z sets land on random signal less often than real or shuffled peptides do.
- That is why every arm, our old base and DIA-NN with its own decoys all sit above nominal. DIA-NN is about 1.2× at
  q ≤ 0.05.
- It is not fixed by scoring features (featD), by interference evidence (decoys share signal as often as
  entrapments), or by adding parents (PAIRLIB).
- Which construction step lowers the random-match rate is not resolved here. Candidates:
  - the precursor mass shift, which moves decoys into other isolation windows with different spectral density;
  - the K/R swap changing the y-ion series;
  - fragment m/z placement relative to common background ions.

## Suggested next step: a frozen plan, which changes the decoy population and so needs approval
- Regenerate the decoys from the same targets so that they match noise like targets do, then rerun PAIRLIB on codon
  (about 1 h). These are not hard decoys: nothing is tuned to resemble particular targets.
- **Arm A:** OpenSwathDecoyGenerator `shuffle`, keeping the termini and the precursor mass. This is the same kind of
  object as an FDRBench entrapment.
- **Arm B:** pseudo-reverse without the K/R swap, which keeps the C-terminal residue and the mass.
- **Label-free pre-check, before any training:** the null-window random-match rates (support ≥ 4, cosine > 0.8) of the
  new decoys must equal the targets' within noise. This needs only the search output.
- **Predictions to freeze:**
  - pooled rho at q ≤ 0.05–0.20 is ≤ 1.15;
  - rule 1c passes on both raws.
