# Voxel v3-seq: peptide sequence encoder (ReDeeM-style atom graph) added to the single decoy-whitened network (written before the featurizer runs and before any training)

User request (2026-11-10): encode the peptide sequence in the same way as the ReDeeM foundation model, and add the encoded features to the voxel model and the heuristic features. ReDeeM predictions and weights are not used.

## Encoder (ported, randomly initialised, trained end to end)

- Featurizer `seq_featurize.py`: ReDeeM residue-local heavy-atom graphs (shared N-CA-C(=O) backbone + side-chain template), 12 atom features (atomic number, mass, electronegativity, valence, formal charge, aromatic, donor, acceptor, backbone, side chain, pseudo-mass, terminal flag), bond-order adjacency with self-loops, exact PTM graphs for UniMod:4 (C) and UniMod:35 (M). These are the only modifications in this library. L = 64, A = 20.
- Encoder (the ReDeeM `PeptideFoundationEncoder` layout, reduced size): atom linear 12→64; 3 graph message layers (row-normalised adjacency, self + neighbour projection, ReLU, output projection, residual + LayerNorm); mean over atoms; linear 64→96; + residue embedding (21) + position embedding (64); LayerNorm; 2 Transformer blocks (4 heads, FF 192, dropout 0.05); LayerNorm; masked mean over residues → 96. Then a linear layer 96→16 and tanh gives the sequence features.

## Arms (the same decoy-whitened objective, schedule, folds and population as PLAN_V3_P2b)

| Arm | Input to g | Head h |
|---|---|---|
| DS11 | [z, sequence] | [z, tanh(g)] |
| DVS | [voxel body, z, sequence] | [z, tanh(g)] |

D11 and DV (PLAN_V3_P2b) are the controls without the sequence input.

## Known risk, measured before this plan (SEQUENCE_TELL.txt, NULL_BASELINE.txt)

- The library sequence alone separates targets from decoys: logistic regression AUC 0.537, human 0.552, entrapment 0.521. The pseudo-reverse decoys keep the C-terminal K/R.
- A model can therefore learn "target-like sequence", which biases the decoy q value optimistically. The entrapment check does **not** detect this fully: the bias favours human over entrapment, which makes FDP_c look better. So a separate, label-free check is needed.
- **Null-window check:** windows with |library RT error| > 3σ (no peptide is expected there; about 178k target and 167k decoy windows per raw). AUC(target > decoy) of the arm score, over null windows and over per-precursor best null windows (the null winners). Baselines (null winners): BASE 0.518 / 0.523, D11 0.502 / 0.502, DV 0.500 / 0.499.

## Pass rule (all on both raws)

1. PLAN_V2 rules 1–4 against the population BASE (ΔN@2 % interval > 0, ΔN@1 % interval not entirely < 0, FDP_c at q ≤ 0.01 ≤ 0.015, entrapment AUC ≤ BASE).
2. **Sequence-tell guard:** null-winner AUC ≤ the control's null-winner AUC + 0.005 (DS11 against D11, DVS against DV). If not, the arm fails, whatever its N.
3. Gain from the sequence input: N@2 % of the arm minus its control > 0 on both raws (reported with the bootstrap interval against the control; a fail of this rule means "the sequence adds nothing").

An arm that passes 1–3 is a candidate for C++ `OpenDIAScoringModel` under a new plan. Otherwise it is reported to the user.

## Run order

DS11 on both raws first (about 10 min). DVS on both raws (about 60–80 min). Evaluation: `evaluate_v2.py` (+ bootstrap of the arm against its control) and `null_auc.py`.
