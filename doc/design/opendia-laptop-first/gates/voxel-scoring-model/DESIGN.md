# OpenDIAScoringModel: one voxel scorer instead of CNN + LDA

Status: design, before any pilot code or run. Raw 03 stays untouched. DIA-NN and OSW are references only.

## Decision

OpenDIA gets **one** scoring model. It replaces the native CNN logit and the stage-2 LDA (with its semi-supervised iterations). It lives in new files, `OpenDIAScoringModel.{h,cpp}` with its own tool, so it cannot be confused with `OpenDIACNN` (CNN, TOPAZ options, stage 2). `OpenDIACNN` stays unchanged as the comparison base until the new model beats it. It is then retired, not merged.

The q-values come from the single model's logit, through paired, pooled target-decoy competition (the existing `-paired_tdc -pooled_tdc` logic, reused unchanged). There is no second model on top.

## Pipeline

1. **Calibration** (`OpenDIACalibrate`, unchanged). The outputs are:
   - anchor calibrants: ID-free targets that beat their paired decoy, with q ≤ 0.01 (513 on raw 01);
   - RT, m/z and mobility calibration;
   - extraction window widths.
2. **Candidate search** (`OpenDIA`, unchanged). Every library assay (target and decoy) is extracted inside its calibrated RT window. Up to 5 candidate peaks are kept per precursor (2.8M windows on raw 01).
3. **Voxel extraction.** One input per candidate window (see "Inputs").
4. **Training set.** Built from calibration and the library only (see "Training"). It never uses search scores, OSW or DIA-NN.
5. **Training and scoring.** Peptide families are split into folds (default 3). A model trains on the training rows of two folds and scores **all** windows of the third. No window is ever scored by a model that saw its family.
6. **Precursor winners and q.** For each precursor, the best window by logit wins, then target-decoy competition gives q.

## Inputs: pilot options a–d

All four share one TOPAZ-derived body:
- a branch per modality: 1D convolutional encoder + co-elution block;
- missing modalities give exact zeros;
- one head that outputs one logit.

Only the input stage changes.

| Option | Input | Change to the body |
|---|---|---|
| a | Today's tensor: 6 fragments × 32 RT (mobility collapsed), MS1 M+0..2 × 32 RT | none (control) |
| b | a + one ppm-offset channel per fragment and RT bin (intensity-weighted, relative to the run's calibrated offset; 0 when there is no signal) | the fragment branch takes 2 channels per fragment |
| c | a with m/z as an axis: 7 sub-bins over ±tolerance for each fragment and RT bin | short convolution over m/z, then the RT encoder |
| d | b with up to 12 fragments and MS1 M−1..M+3, with a validity mask | channel count becomes a setting, with masked co-elution |

Mobility, when present (timsTOF), is a further axis on the same pattern. On the Astral design raws it is absent and contributes zeros.

## Training

| Group | Source | Label | Weight |
|---|---|---|---|
| Positives | Calibration anchors, at the candidate window that contains the anchor's run RT | 1 | 1 |
| Paired decoys | Each anchor's paired decoy, at its best candidate window | 0 | 1 |
| Random decoys | Decoy candidate windows from all over the library | 0 | 1 |
| **Decoy-like negatives** | Target candidate windows that are very unlikely to be the peptide (rules below) | 0 | 0.5 |
| **Synthetic false targets** | Built from a decoy or decoy-like window by adding a scaled fragment trace of an unrelated present window to some of its channels (partial, shifted or wrong-ratio co-elution) | 0 | 1 |
| Synthetic weak positives | An anchor window scaled down (s in [0.05, 0.5]) and added onto a decoy-like window | 1 | 1 |

**Decoy-like negatives.** The aim is target windows whose signal looks like evidence but is not the peptide. That is the distribution that contains the false targets. A target window qualifies if both of these hold:
- it is a non-winning candidate of its precursor (that precursor has a higher-scoring window by the search's own transition support, which is label-free); or its centre lies more than 3σ of the calibrated RT residual from the predicted RT (a calibration quantity, not a score);
- its precursor is not an anchor.

About 0.1 % of random target windows are true peaks (estimated earlier). The RT rule removes most of those. The weight of 0.5 limits the effect of the rest.

**Why these negatives matter.** In rounds 1–6 the only negatives were decoys. A deep model then learns "real peptide-like signal against decoy-like signal". That lets entrapment winners (real sequences matched to interference) rank above decoys (entrapment-vs-decoy AUC 0.52). Decoy-like target negatives and synthetic false targets carry the same "peptide-like signal" features as false targets, so a model can only separate them by evidence that **this** precursor is present: mass accuracy, co-elution of all fragments, library ratios and MS1 isotopes.

## Leakage rules

- Entrapment status is never read by any training script. It is used only in the evaluation file.
- Anchor selection uses calibration only. Decoy-like selection uses the search's label-free statistics and calibration only.
- Fold assignment is by peptide family (unmodified parent sequence), so a target and its decoy are always in the same fold.

## Pilot (offline, Python, raw 01 only)

- **Data.** `OpenDIA -out_xipm` is used for the training groups and the evaluation sample (`astral-voxel-pilot-20261110/xipm_precursors.txt`). The evaluation sample has base q ≤ 0.01 targets, entrapment winners, decoy winners and weak targets.
- **Grid.** Inputs a–d × two training sets:
  - T0: anchors + decoys only;
  - T1: the full table above.
- **Measured on held-out folds** (reported per cell):
  - AUC of present over decoy;
  - AUC of entrapment over decoy;
  - N at FDP_c 1 % and 2 % on the evaluation sample (a proxy; the sample is not the full run).
- **Pass**, for a pilot cell to move to C++:
  - entrapment-over-decoy AUC within 0.5 ± 0.01, AND
  - present-over-decoy AUC ≥ that of cell (a, T0), AND
  - N at FDP_c 2 % above (a, T0), with a bootstrap 90 % interval of the difference above 0.

This rule is fixed before any pilot training. It is copied, with a hash, into `PLAN_PILOT.md` before the first model is trained.

## C++ (after the pilot)

`OpenDIAScoringModel` is a new library file and a new tool, built on the libtorch helpers that `OpenDIACNN` already uses. Defaults elsewhere stay byte-identical (HeLa check). The full-run gate on raws 01 and 02 uses the existing `PLAN.md` rules 1–4, against base CNN + LDA 1029 / 1298. Raw 03 needs approval.
