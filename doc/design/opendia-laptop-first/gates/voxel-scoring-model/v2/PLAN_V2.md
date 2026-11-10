# PLAN_V2: voxel model with heuristic features and semi-supervised rounds

Written before any v2 training. The user asked for two changes on 10 November 2026, after the voxel pilot failed (`../RESULT_PILOT.md`):
1. Add the heuristic features to the voxel model, including the b/y-ion counts.
2. Semi-supervised rounds for the voxel model.

A lower weight for decoy-like negatives comes later and is not in this plan.

## Is the model semi-supervised now?

- The **base** is semi-supervised in its second stage only. `OpenDIACNN -stage2` refits the LDA in 3 rounds per fold, with positives at q ≤ 0.15, then 0.05. The CNN itself trains once on calibration seeds (`-ss_iterations 0`).
- The **voxel pilot** was not semi-supervised. It trained once on about 500 anchors and their decoys.
- The **single TOPAZ rounds 1–5** used `-ss_iterations 2`. That made the result worse there (149 at q ≤ 0.01 after round 1, CALIB §13). That model had no early-stopping fold, and it ran on 600 seeds. This plan uses fresh models per round, a validation fold, and the larger positive sets of the training folds.

## Data

- **Design raws 01 and 02.** Raw 03 is not used.
- **Population** (`make_population.py`, label-free): every candidate window of every target-decoy pair whose best `transition_support` is ≥ 3.
  - raw 01: 94,109 pairs, 186,702 precursors, 1,270,510 windows;
  - raw 02: 95,447 pairs, 189,464 precursors, 1,275,383 windows.
  - Reason: the raw-point export and the voxels of all 2.8M windows do not fit on this laptop.
  - Effect on the base, measured before training, with the same TDC and metric as the arms: N@FDP_c 1 % is 1024 (raw 01) and 1288 (raw 02), against 1029 and 1298 on the full run. This **population base** is the comparison for every arm.
- **Inputs:** voxel option b from the pilot (fragment XIC, ppm channel, MS1 M+0..2, library intensities), from `OpenDIA -out_xipm`. The export's candidate table must be byte-identical to the base search table (`export0N/CMP.txt`).
- **Heuristic features:** the 11 base stage-2 features (search 4, boundary 5, byion 2), read from the base `features.tsv`, standardised with the training folds' mean and SD.
- **Folds:** 3 folds by peptide family. Each window is scored only by the model of the fold that did not train on its family.

## Arms (`train_v2.py`)

| Arm | Heuristic features | Semi-supervised rounds |
|---|---|---|
| V0 | no | no (control, the pilot cell b_T0 on the full population) |
| F0 | yes | no |
| S0 | no | 2 |
| FS | yes | 2 |

- **Round 0:** calibration anchors, their paired decoys and 6,000 random decoys (pilot T0).
- **Rounds 1 and 2:** score the training folds, take the pooled paired TDC winners, and train a new model on target winners at q ≤ 0.15 (round 1), then 0.05 (round 2), against decoy winners of the training folds (4× the positives, at least 6,000).
- Training is otherwise as in the pilot: Adam 1e-3, batch 64, at most 40 epochs, early stopping with patience 5 on 15 % of the families.
- Entrapment labels are read only by `evaluate_v2.py`.

## Decision rule (per arm, against the population base, on both raws)

The metric is N at entrapment FDP_c 1 % and 2 % (FDP_c = N_e (1 + 1/r) / N), on precursor target winners after pooled paired TDC. The paired bootstrap uses 1,000 precursor resamples, seed 5. This rule uses the bootstrap because `metric_noise.txt` showed that a fixed N threshold is noise-limited. **This is a rule for this offline plan only.** The gate for the C++ tool stays `PLAN.md` rules 1–4, unless the user approves a change.

An arm **passes** when all of these hold on **both** raws:
1. The bootstrap 90 % interval of N@2 % minus the population base lies above 0.
2. The bootstrap 90 % interval of N@1 % minus the population base does not lie entirely below 0.
3. FDP_c at q ≤ 0.01 is ≤ 0.015.
4. The entrapment-over-decoy winner AUC is ≤ the population base AUC (0.531 on raw 01, 0.535 on raw 02).

An arm that passes goes to C++ `OpenDIAScoringModel` and then to the full gate under a new plan. If no arm passes, the result is recorded, and the next step is the user's decision.

**Reported only:** q ≤ 0.05 N and FDP_c, N at FDP_c 0.5 % and 1.5 %, positives per round, epochs, time and memory.
