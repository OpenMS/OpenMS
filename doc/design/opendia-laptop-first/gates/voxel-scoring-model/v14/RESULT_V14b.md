# RESULT_V14b: Local-flank weak-profile counterfactual

**Run:** raw01 only, 12 November 2026  
**Plan:** `PLAN_V14b.md`, frozen at 2026-11-12T00:40:00Z; SHA256 `506b3b852f0b6eab8cc9f6594eab6234e41dff5a901de314bcaf3547ee8a9caa`  
**Code:** `weak_aug_v14_features.py` SHA256 `4b9cee8d0a0f24a39c3a98c3eeae45d2e665de3b6f1c58af9db3943486370119`; instrumented `train_v4.py` SHA256 `a3a18aa5598a715a3b8aadd075db6a800e98143cd4fb255c03f39a520c2c72a4`; evaluator SHA256 `0cc8e77fd48aa4964b55c55c0e68b1b7706f4feaf8f854374362e3de4bf79446`.

## Question and method

Does D11FWT reduce the score of an otherwise accepted, coherent peak group when its positive XIC signal above the measured local flank baseline is attenuated to the observed Stage-D/accepted DIA-NN area ratio? The factor was fixed at 0.25655589613105056 from the previously recorded raw01 median log-area difference. The baseline was each fragment trace's median in bins 0:8 and 24:32. The existing D11F model was reconstructed with its unchanged training recipe, and the counterfactual six F12 summaries were scored by the held-out fold scorer. No intervention model was trained.

This is deliberately a partial-input diagnostic: the 11 heuristic features were left at their observed values. Only the six F12 inputs were recomputed. The feature-generation code loaded no entrapment, DIA-NN, or OSW labels. Evaluation labels were joined only after counterfactual scores had been written.

## Integrity checks

- 1,365 raw01 target-decoy winners passed the existing pooled q<=0.005 donor gate. All had complete XIC profiles.
- Recomputed F12 values exactly matched the stored values on the validation subset (maximum absolute difference 0).
- Reconstructed D11F scores had Spearman rho 1.0 against the saved D11F scores.
- W-selected row identity agreement was 1.0.
- Reapplying the saved per-fold T calibration reproduced finite D11FWT scores with maximum absolute error 1.91e-6.

## Result

| Measure | Result |
|---|---:|
| Baseline donor median D11FWT score | 31.4731 |
| Attenuated donor median score | 30.8034 |
| Median paired calibrated score change | -0.000004 |
| Stage-D median score | 5.6664 |
| Stage-E median score | 30.8690 |
| Attenuated donors inside observed Stage-D score range | 2.64% (baseline 0%) |
| Attenuated donors below Stage-E 10th percentile | 7.77% |
| Target-vs-same-pair-decoy AUC, counterfactual donors | 0.99997 |
| Baseline entrapment-vs-decoy AUC | 0.49067 |

The median synthetic donor remains at the Stage-E score level, not the Stage-D level. Attenuation moves only 2.64% of this high-confidence donor set into the observed Stage-D score range. The large same-pair AUC is selection-conditioned: donors were selected from baseline winners at q<=0.005, so it is not evidence that the transformation recovers Stage-D targets or would preserve global error control. The entrapment control AUC remains below 0.53, but it is unchanged baseline evidence, not a score for a trained augmentation.

The six F12 summaries are mostly invariant to positive signal attenuation: median changes are 0 for support, co-elution, library cosine, extra support, and extra co-elution. The only summary with a nonzero median change is `f12_extra_share` (+9.31e-5). Although every row has some numerical change, the change is too small for a median score movement. This is consistent with F12 mostly describing which fragments are present and whether their profiles co-vary, rather than absolute peak amplitude.

## Conclusion

**No evidence that weak-profile augmentation through the current F12 input path explains or recovers the Stage-D misses.** The existing network cannot infer a substantial amplitude reduction from features that largely discard amplitude. This does not prove that the whole D11FWT model is amplitude-invariant, because its 11 heuristic features were not regenerated from attenuated profiles. It also does not show that adding synthetic weak examples to training would help. Do not train an augmentation model, change the base/pass rule, or promote C++ from this result.

The paired-decoy V14 attempt was stopped before feature extraction: q<=0.001 produced no donors because the add-one q minimum was about 0.001499, and only 531/1,327 paired decoy windows were within 5 seconds at q<=0.005. V14b replaced that infeasible matching assumption with each trace's local flank baseline. The original V14 plan and hash remain untouched.

Artifacts are in `v4/weak_aug_v14/`: counterfactual F12 NPZ, fold score NPZ files, `counterfactual_scores_raw01.tsv`, `FEATURES_STATUS.txt`, `EVALUATION.txt`, and execution logs. The unchanged D11FWT reference arrays were not overwritten.
