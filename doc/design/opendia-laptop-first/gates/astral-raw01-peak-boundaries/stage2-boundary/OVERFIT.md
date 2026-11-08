# Overfitting and error-control checks for a single CNN → MLP discriminator (raw 01, design only)

This is a post-hoc diagnostic on raw 01, a reused design run. It is not a gate.

## What was measured

Entrapment labels come from the library protein accessions only (`_p_target`), with r = 0.9851, and were never used for fitting. FDP_c is the "combined" estimate from Wen et al. 2025. All counts are paired-pooled q, offline.

| Discriminator (search+boundary features unless noted) | q ≤ 0.01 N | N_e | FDP_c | q ≤ 0.05 N | FDP_c |
|---|---|---|---|---|---|
| Tool stage 2 (LDA) | 746 | 3 | 0.008 | 1441 | 0.070 |
| LDA, Python re-implementation | 745 | 3 | 0.008 | 1440 | 0.070 |
| MLP 32-16, 5-seed ensemble | 1009 | 9 | **0.018** | 1534 | 0.070 |
| MLP 16 + dropout 0.2 + weight decay 1e-2, 5-seed ensemble | 935 | 6 | 0.013 | 1524 | 0.093 |
| MLP, nested out-of-fold label selection, 5-seed ensemble | 1111 | 10 | **0.018** | 1602 | 0.073 |
| MLP on CNN + 4 search scalars only, 5-seed ensemble | 654 | 2 | 0.006 | 1526 | 0.066 |

Other checks:

- **Fold hygiene:** 0 of 287010 target/decoy pairs span folds.
- **Memorisation probe:** this is the AUC of held-out decoys against the decoys a fold model trained on, scored by that model, where 0.5 means no memorisation.
  - LDA: 0.50 / 0.50 / 0.51.
  - MLP: 0.60 / 0.55 / 0.57.
  - Regularised MLP: 0.63 / 0.58 / 0.65.
- **Poisson check:** if the true FDP were exactly 1%, seeing N_e ≥ 9 among 1009 has probability 0.07, and N_e ≥ 10 among 1111 has probability 0.05.

## Reading

1. **Cross-fitting already protects against classic train-on-test overfitting.** The MLP does memorise its training decoys, but scoring is out-of-fold and the pairs never leak.
2. **The real risk is different.** A flexible model can learn what a decoy looks like rather than what a wrong match looks like. When it does, decoys are pushed down but false targets are not, so the decoy count says 1% while entrapment says about 1.8%. The R symmetry check cannot see this, because it uses decoys only. On raw 01 every MLP variant sits at or above the 1.5× pass line used by the entrapment gate. LDA sits at 0.8%. With N_e between 6 and 10, the evidence for a real excess is borderline (p ≈ 0.05–0.07). It is not proven.
3. **More regularisation was not the fix.** One predeclared variant was tried, with no sweep. It memorised more, kept less gain, and was worse at 5%.
4. **The gain does not come from the CNN.** The CNN plus search scalars through the MLP gives 654 at 1%. The gain needs the peak-group scalars (co-elution, library correlation, apex agreement) combined non-linearly.

## Implication for the design

One model is the right shape: CNN trace embedding plus a fixed set of always-computed scalars, going into a small MLP head, trained cross-fitted. This is the TOPAZ `CandidateScorer` design, and it removes the separate stage-2 and opt-in feature sets.

It must ship with:

- error control that does not rely on decoys alone (an entrapment FDP check on an unused run, as the acceptance rule);
- pair-grouped folds;
- out-of-fold label selection;
- a seed ensemble;
- a predeclared pass rule that includes FDP_c ≤ 0.015.

Candidate mitigations for decoy-artefact learning, to be decided before the run:

- monotone or sign-constrained inputs for scalars with a known direction;
- excluding inputs whose decoy distribution is known to differ from false targets (decoy-generation artefacts);
- training against decoys only on the non-artefact inputs.
