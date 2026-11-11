# Voxel v12: reference-profile interference correction as a second view for the learned voxel model

## Purpose

Test whether a fragment-group reference profile and interference-corrected chromatogram view can help recover weak stage-D precursors without letting the correction manufacture confidence. This is a new, bounded model experiment. It does not change the base or pass rule.

## Evidence and pre-plan look

- **Base:** D11FWT, old population, same raw01/raw02 inputs as V9-V11.
- Previous raw-spectrum scans counted fragment matches but did not model trace-level interference. The stored arrays have six fragment XICs over 32 RT bins. The existing D11FWT run in `train_v4.py` sets `voxel=False`; it does not learn directly from these arrays. It uses engineered features, including the 11 search heuristics and the six-fragment-derived feature group F.
- **Raw01 pre-plan look only** (`INTERFERENCE_DIAGNOSTIC_raw01.txt`; score-matched using D11FWT strata):
  - D-vs-decoy AUC / entrapment-vs-decoy AUC: reference quality 0.584 / 0.472; raw profile correlation 0.597 / 0.467; corrected library cosine 0.641 / 0.541; corrected core fraction 0.642 / 0.522.
  - All-human-vs-decoy AUC / entrapment-vs-decoy AUC: corrected library cosine 0.584 / 0.549; corrected core fraction 0.553 / 0.540.
  - These are modest candidate signals, not evidence of a likely win. Entrapment is not strongly enriched by these features, but is slightly above 0.5 for some features. The full model and the frozen entrapment guard remain decisive.
- **Not looked at:** raw02 corrected profiles or any trained v12 score.

## Mechanism

For each candidate window, the label-free preprocessing computes a fragment reference trace as the active trace with the highest average Pearson correlation to the other active traces. It smooths sparse binned traces with a centered three-bin moving average. Other traces are scaled to the reference using reference-squared weighted intensity; for poorly correlated traces (r < 0.8), the scale uses the minimum trace/reference ratio near the reference apex. Each smoothed trace is capped at 1.5 times its scaled reference. This is DIA-NN-inspired, adapted to these sparse 32-bin arrays, and is not claimed to reproduce DIA-NN's implementation.

**Do not replace the raw profiles.** Save both the raw and corrected profiles. Do not compute or train on corrected-profile coherence, because the cap itself forces the corrected trace toward the reference shape. Feed both profile views to one learned voxel model, alongside the existing D11FWT inputs; the network can learn whether, where, and how much to use the correction. Correction strength and corrected-view summaries are retained for diagnosis, not appended to the model in this plan.

## Arm

**D11FIWT:** existing D11FWT inputs and training procedure, plus the dual-view voxel body consuming raw six-fragment traces, corrected six-fragment traces, ppm maps, MS1 profiles, and library intensities. `I` enables the voxel body; `F` retains the existing six-fragment feature group. The corrected view is generated without labels. The model remains one cross-fitted scorer with the existing within-precursor W step and per-fold decoy-tail transform.

## Frozen pass rule

Both raws must satisfy the same v10/v11 rules:
1. PLAN_V2 rules 1-4 against BASE, including FDP_c at q <= 0.01 <= 0.015.
2. N@2% greater than D11FWT on both raws.
3. Null-window winner AUC no more than D11FWT + 0.005.

Report N at 0.5/1/1.5/2%, exact q<=0.01 FDP_c, entrapment AUC, null-winner AUC, stage funnel, and paired bootstrap versus D11FWT. An arm that passes becomes the base. Otherwise retain D11FWT and stop this correction approach.

## Implementation and validation

- `interference_profiles.py` writes `xic_corrected.npy` and `interference_features.npy`; original `xic.npy` remains unchanged.
- `train_v4.py` enables the two-view voxel path only when `I` is present in the arm name. Legacy arms retain their previous feature-only path.
- Synthetic checks cover an added off-peak interference bump, all-zero traces, and a single active trace. Model forward-shape check verifies a 200-value body embedding.
- Run `chain_v12.sh`, which checks this plan's sha before generating raw02 features or training.

## Risks and limits

- The reference can itself be interfered with; all fragments can share a true peak shape; the correction may remove legitimate signal. Keeping the raw view is intended to let the scorer reject harmful corrections.
- The summary features are correlated with existing coelution and library-fit measures. The primary novelty is the raw-plus-corrected time-series input, not the seven diagnostic summaries.
- Both raws are the only validation set for this work. A pass is not general validation beyond ProteoBench Astral raws 01 and 02.
