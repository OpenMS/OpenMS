# Stage-2 feasibility on raw 01 (design run, target/decoy labels only)

This follows the raw 02 post-mortem. The question is whether extraction-level evidence (fragment co-elution, library intensity correlation, MS1 mass error) adds target-decoy separation beyond the CNN logit. It uses raw 01, the design run, and no entrapment labels, per `../PROTOCOL.md`. It is a design check, not a gate.

## Steps

1. `select.py`: from the frozen design scores (`../pipeline/crossfit_it0`, iteration 0, decoy-standardised as in `-pooled_tdc`), take up to 60 target and 60 decoy precursor winners per CNN z band.
   - Only 8 decoys exist above z 10, so that band has 8.
   - A first pass by mistake used `crossfit_pooled` (`-ss_iterations 1`, the overfit run). It was discarded and kept in `wrong_ss1/`.
2. `run.sh`: the frozen OpenDIA search with `-out_xipm` for these 548 precursors.
   - The search output is byte-identical to `../pipeline/search.tsv`.
   - 4932 peak maps and 556939 points, in 2.5 min and 4.4 GB.
3. `features.py`: per precursor, within ±4 s of its CNN winner apex, computes:
   - share of fragments with signal;
   - mean pairwise fragment-trace correlation;
   - sqrt-intensity correlation with library intensities;
   - fragment and MS1 ppm error;
   - MS1-to-fragment co-elution.

## Result (`feasibility.txt`)

Medians, target / decoy:

| CNN z | co-elution | library corr | abs MS1 ppm |
|---|---|---|---|
| ≥ 10 | 0.89 / 0.90 | 0.93 / 0.74 | 1.3 / 7.1 |
| 8 to 10 | 0.74 / 0.28 | 0.85 / 0.52 | 2.1 / 7.3 |
| 6 to 8 | 0.32 / 0.17 | 0.63 / 0.40 | 4.7 / 7.7 |
| 4 to 6 | 0.09 / 0.05 | 0.17 / 0.21 | 7.2 / 6.4 |

- **Where the information is:** in CNN z 6 to 10, which lies just below the current 1% cut, co-elution, library correlation and MS1 mass error all still separate targets from decoys.
  - MS1 ppm error is the strongest single feature. It also separates the 8 top decoys, which co-elute well (0.90) but sit 7 ppm off in MS1.
- **Hand-chosen rule** (co-elution > 0.6, library correlation > 0.6, |MS1 ppm| < 4): the sampled pass rates times the band populations give about 690 targets against about 3 decoys above z 6.
  - The current CNN-only result on this run is 321 at pooled q ≤ 0.01.
  - This extrapolation rests on 60 samples per band, and the rule was picked by looking at raw 02. It is therefore only a feasibility signal, not a design result or an ID count.

## Decision

A native, opt-in second-stage score is worth building. It would:
- compute these features for every CNN window from the extraction OpenDIACNN already performs;
- combine them with the CNN logit through a semi-supervised linear discriminant, PyProphet-style, trained on cross-fit folds;
- then run the same pooled paired TDC.

Design will be on raw 01 only, using target/decoy counts. Evaluation needs raw 03, the last held-back replicate, under a new hashed gate.
