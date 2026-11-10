# Voxel v3 phase 2b: one network with a decoy-whitened (deep LDA) objective (written before any run)

Motivation (RESULT_P2a.md): objectives that use the positive-class spread (BCE, the Fisher loss) prefer entrapment. Positives are a mix of true and false targets. The decoy-only covariance read-out kept the IDs and had the lowest entrapment AUC (post-hoc). This phase makes that read-out the objective of **one** network. There is no second-stage model.

## Model (one network, one score)

- Input: z = the 11 heuristic features, standardised on the training folds. Arm DV adds the voxel tensors (pilot option b: fragment XIC + ppm channel, MS1 M+0..2, library).
- Representation: h = [z, tanh(g(u))]. g is a 2-layer MLP (64 hidden, 16 out). u = z (D11) or [voxel encoder parts (v2 Scorer body, 176 values), z] (DV). Arm D11n uses h = tanh(g(z)) only (no z skip).
- Score: s = w·h, w = (S₋ + λI)⁻¹(μ₊ − μ₋), with S₋ the covariance of the negatives only and λ = 1e-2.
- Loss on each batch (256 positives + 1024 negatives, w solved in closed form on the batch): −(μ₊ − μ₋)ᵀ(S₋ + λI)⁻¹(μ₊ − μ₋).
- Final w for scoring: solved once on all training rows of the round.
- Adam 1e-3, max 30 epochs (100 steps per epoch for D11 / D11n, 50 for DV), early stopping on the validation-family objective (patience 5).

## Schedule, population, metrics

As phase 1 / L11d: 3 family folds, round 0 on T0 (anchors against anchor decoys + random decoys), rounds 1–2 at q ≤ 0.15 / 0.05 (target winners against all decoy winners of the training folds), a fresh model per round, family-split validation. Population, evaluation and bootstrap: `evaluate_v2.py`. Entrapment labels are not used in training.

## Gate 1: 11 features (D11; D11n is for information only)

D11 passes gate 1 if, on both raws: N@1 % ≥ 0.95 × L11d (1259 / 1244) and the entrapment AUC ≤ L11d + 0.005 (0.5274 / 0.5273). Note: λ = 1e-2 here, and L11d used 1e-3.

## Gate 2: voxel + 11 features (DV), runs only if D11 passes gate 1

DV passes if all PLAN_V2 rules hold on both raws against the population BASE: (1) the ΔN@2 % 90 % interval is above 0; (2) the ΔN@1 % interval is not entirely below 0; (3) FDP_c at q ≤ 0.01 ≤ 0.015; (4) the entrapment AUC ≤ the BASE AUC. DV − D11 is reported to show what the voxel input adds.

A DV pass is a candidate for C++ `OpenDIAScoringModel` under a new frozen plan with the full gate. A fail is reported to the user.
