# Voxel v3 phase 2a: an objective for one network that does not prefer entrapment (written before the run)

From RESULT_P1.md: with the 11 heuristic features, BCE-trained models (linear or MLP) prefer entrapment (AUC 0.554–0.558). The LDA does not (0.531). This phase tests objectives for a single network. Same population, folds, T0 rows, rounds (q 0.15 / 0.05) and evaluation as phase 1. Inputs: the 11 features only.

| Arm | Head | Objective | Negatives in rounds |
|---|---|---|---|
| G11a | linear | BCE (class-balanced) | all decoy winners (separates the sampling effect from G11) |
| FL11 | linear | Fisher loss: −(μ₊ − μ₋)² / (σ²₊ + σ²���), computed on balanced batches (1024 positives + 1024 negatives) | all decoy winners |
| FM11 | MLP (as N11) | Fisher loss as FL11 | all decoy winners |

Training for Fisher arms: Adam 1e-3, 300 steps per epoch, max 40 epochs, early stopping on the validation Fisher loss (patience 5).

Reading:
- An arm is **acceptable** if, on both raws, N@1 % ≥ 0.95 × L11 and entrapment AUC ≤ L11 AUC + 0.005.
- FL11 must be acceptable: it is the gradient form of the LDA. If not, the implementation is wrong; stop and report.
- FM11 acceptable → phase 2b puts the voxel encoder under the Fisher objective (one network, voxel + features).
- G11a acceptable → sampling was the cause; phase 2b uses BCE with all negatives.
- No flexible arm acceptable → report to the user.
