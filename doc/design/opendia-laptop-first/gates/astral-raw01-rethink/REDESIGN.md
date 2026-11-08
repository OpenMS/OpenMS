# Rethink: what limits OpenDIA on Astral raw 01, and a redesign (design diagnostic)

Raw 01 is a design run, and the user approved reusing raws 01 and 02. DIA-NN 2.3.1 is a reference only. Its IDs are used here to measure the ceiling and to simulate a predictor. They are never used as a pipeline signal.

## Ceiling and funnel (`ceiling01.py`)

The reference set is the DIA-NN 1% precursors that exist exactly in our library: 2754 on raw 01, 11 of them entrapment.

| step | precursors |
|---|---|
| any search candidate | 2741 |
| a candidate within 10 s of DIA-NN RT | 2449 |
| our winner within 10 s | 2290 |
| accepted at q ≤ 0.01 (tool, search+boundary) | **730** |
| accepted at q ≤ 0.05 | 1297 |

Recall at 1% by intensity tertile is 1% / 23% / 71% (low / mid / high).

## The current features are the main limit, not only the learner

The oracle here is fully supervised (DIA-NN-confirmed right peaks against decoy winners). It is cross-fitted, uses all 17 current features, and gives an upper bound for any learner on these inputs.

| oracle | decoy/target ratio 1% | 5% |
|---|---|---|
| LDA | 893 | 1336 |
| MLP 32-16 | 1091 | 1580 |

Even with perfect labels, the current features recover about 1100 of the 2449 findable precursors at 1%. A better combiner can move us from 730 towards about 1100, but not towards 2449. The missing part is evidence the features do not carry, mostly for the low and mid intensity tertiles.

## RT is the largest missing piece of evidence (`rtlever_full.py`)

- Our library RT, after the best non-linear global map to DIA-NN RT, has a robust σ of **24 s** (q90 45 s). DIA-NN's in-run predicted RT has σ **0.6 s**. Peaks are about 10 s wide, so our RT prior fits about 5 peak widths into ±1σ.
- **Simulated lever:** the oracle score plus a Gaussian RT term, with candidates limited to ±3σ, over the full 287K-pair TDC. The DIA-NN-confirmed precursors get a simulated predictor of precision σ, and everything else keeps the calibrated library RT.

  | RT σ | accepted at q ≤ 0.01 (confirmed, right RT) | q ≤ 0.05 |
  |---|---|---|
  | none | 1134 (1088) | 1873 |
  | 24 s (current) | 1219 (1157) | 1846 |
  | 12 s | 1358 (1304) | 2028 |
  | 6 s | 1527 (1475) | 2078 |
  | 3 s | 1627 (1587) | 2153 |

  This is optimistic, because only the confirmed set gets the better RT. It does show that RT precision is worth more than any combiner change: +37% from 24 s to 3 s on top of the oracle.
- **Teacher-free refit, first try** (`rtrefit.py`): our 746 accepted IDs, ridge on sequence composition and a calibrated-RT cubic, cross-fitted. σ went from 29.4 to 25.7 s on held-out confirmed precursors. A composition model is too weak. The remaining error is per-peptide error in the library prediction, and a global map cannot remove it.
- **Pretrained redeem `rt_cnn_tf`** on the same precursors gives σ 137 s against 25.6 s for the library RT, with correlation 0.78 to the library RT. This is probably an input-format or normalisation mismatch and has not been investigated. Fine-tuning on our 614 IDs crashed in redeem at epoch 0: `save_epoch_checkpoint` computes `epoch - 1`, which underflows in a debug build.

## Redesign

The redesign makes one pipeline, iterative like DIA-NN, with no opt-in stages:

1. **Search** (unchanged core) produces candidates and trace windows.
2. **One scorer:** CNN trace encoder (fragments and MS1) plus a fixed set of always-on scalars (search, boundary, extraction and RT likelihood), with a small MLP head (the TOPAZ `CandidateScorer` shape). It is cross-fitted with pair-grouped folds, out-of-fold label selection and a seed ensemble. This replaces stage 2 and the `-stage2_feature_set` options.
3. **RT refinement loop:**
   - take the confident IDs from pass k;
   - fine-tune a sequence RT predictor (redeem `rt_cnn_tf`, in-run, teacher-free, cross-fitted so a precursor never predicts itself) and refit the global map;
   - recompute the RT term and window;
   - rescore.

   Two or three passes. This is the lever the simulation says is the largest.
4. **More evidence for weak precursors** (follow-up): more fragments per precursor (the library has 6; regenerate with 12 and predicted intensities) and MS1 isotope traces in the CNN.
5. **Error control as an acceptance rule:** paired TDC, plus the entrapment FDP_c ≤ 0.015 at 1% on the evaluation run, plus seed variance. `OVERFIT.md` showed decoy-only checks miss a 1.8% FDP_c.

**Target:** about 1500 or more confirmed at 1% on raw 01 and 02 (the OSW+Percolator same-library level is about 1400 confirmed, at FDP_c 3.3%), with FDP_c ≤ 1.5%.

## Experiment order (each one predeclared, design on raws 01 and 02, raw 03 for the gate)

1. **RT predictor sanity:** fix the redeem input or normalisation and the epoch-0 crash, then measure pretrained and fine-tuned σ on held-out confirmed precursors. Pass: σ ≤ 12 s.
2. **Offline end-to-end:** current features plus the RT term from the fine-tuned predictor, one MLP scorer, entrapment FDP_c on raws 01 and 02.
3. **Native implementation** of the unified scorer and the RT loop, removing the stage-2 opt-ins.
4. **Fragments and MS1 isotopes**, if steps 1–3 plateau below target.
