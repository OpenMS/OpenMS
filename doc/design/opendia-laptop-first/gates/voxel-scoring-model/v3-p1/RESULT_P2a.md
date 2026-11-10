# Voxel v3 phase 2a result: no single-network objective is acceptable. The cause is found.

Plan: `PLAN_V3_P2a.md` (sha `96d4a772…`, written 14:59:41Z, checked after the run). Deviations: `PLAN_V3_P2a.deviations` (D1 signed Fisher loss, D2 FM11 stopped with zero positives and counts as a failure). Inputs: the 11 heuristic features only.

| Raw | Arm | N@1 % | N@2 % | FDP_c q ≤ 0.01 | Entrapment AUC | Acceptable (≥ 0.95 × L11 and AUC ≤ L11 + 0.005) |
|---|---|---|---|---|---|---|
| 01 | L11 (LDA, reference) | 1267 | 1484 | 0.010 | 0.5310 | |
| 01 | G11 (linear BCE, sampled negatives) | 1300 | 1345 | 0.004 | 0.5577 | no (AUC) |
| 01 | G11a (linear BCE, all negatives) | 1298 | 1468 | 0.005 | 0.5482 | no (AUC) |
| 01 | FL11 (linear, Fisher loss) | 1118 | 1240 | 0.008 | 0.5529 | no |
| 01 | FM11 (MLP, Fisher loss) | stopped, 0 positives in round 1 | | | | no |
| 02 | L11 | 1261 | 1380 | 0.008 | 0.5306 | |
| 02 | G11 | 1146 | 1242 | 0.002 | 0.5571 | no |
| 02 | G11a | 1218 | 1288 | 0.011 | 0.5487 | no (AUC) |
| 02 | FL11 | 907 | 1173 | 0.009 | 0.5530 | no |

By the plan, FL11 is not acceptable, so we stop and report. The plan assumed "FL11 = LDA". A closed-form diagnostic after the run (`EVALUATION_P2a_diag.txt`, not a frozen test) shows that this assumption was wrong. The FL11 implementation is correct.

## Diagnostic: the covariance weighting is the cause

Same closed-form ridge LDA and same rounds. Only the covariance in w = S⁻¹(μ₊ − μ₋) changes:

| Raw | S | N@1 % | N@2 % | Entrapment AUC |
|---|---|---|---|---|
| 01 | pooled (S₊ + S₋)/(n − 2), the base LDA. About 98 % decoys | 1267 | 1484 | 0.5310 |
| 01 | balanced S₊/n₊ + S₋/n₋ (= the Fisher loss) | 1113 | 1259 | 0.5527 |
| 01 | decoys only S₋/n₋ | 1259 | 1496 | **0.5274** |
| 02 | pooled | 1261 | 1380 | 0.5306 |
| 02 | balanced | 917 | 1136 | 0.5524 |
| 02 | decoys only | 1244 | 1432 | **0.5273** |

- The balanced LDA reproduces FL11 almost exactly (1113 / 917 against 1118 / 907, AUC 0.5527 / 0.5524 against 0.5529 / 0.5530). So the Fisher-loss code is correct.
- L11 was re-run and gives byte-identical scores (cmp).

## Mechanism (well supported on the linear models)

The positive training class is a **mixture**: true targets plus false targets (the round labels are target winners at q ≤ 0.15 / 0.05). The direction that separates true from false targets has high variance **inside the positive class**. Every objective that uses the positive-class spread (S₊ in the Fisher loss, the boundary positives in BCE) gives this direction a low weight, or it fits the false positives. What stays is the "looks like a target entry, against the decoys" direction, and entrapment entries share it. The decoy-whitened mean shift S₋⁻¹(μ₊ − μ₋) measures distance from the null (decoy) distribution only. The false positives then shrink the mean shift, but they do not turn its direction. This is what target-decoy scoring needs.

## Consequence for one voxel network (next phase, needs a frozen plan)

The objective must be **decoy-whitened and must not fit the positive class**. Candidate, one network, no second stage:
- the network produces a feature vector h (voxel encoder + the 11 features);
- the score is a linear read-out w·h, with w = (S₋ + λI)⁻¹(μ₊ − μ₋) solved in closed form on each batch (a differentiable "deep LDA" layer with decoy-only covariance);
- the loss is the decoy-whitened mean shift −Δμᵀ S₋⁻¹ Δμ, with a bound on h (layer norm) so that a few positives cannot inflate the mean.

Test order: on the 11 features first (it must reproduce L11d: N ≥ 0.95 × L11, AUC ≤ L11 + 0.005), then with the voxel encoder.
