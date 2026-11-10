# OpenDIA rewrite: step board

[Master plan](README.md)

**Status, 10 November 2026.** This file lists the steps of each work package and the state of each step. Each row links to the record of the step. The script [`board_sync.py`](board_sync.py) copies the R-WP sections of this file to the GitHub project [OpenMS/38](https://github.com/orgs/OpenMS/projects/38). Edit this file first, then run the script. The "Now" table repeats rows from the R-WP sections, and the script does not copy it.

**States:**
- **Done (pass)**: the step met its predeclared rule.
- **Done (fail)**: the step ran under a predeclared rule and failed. The record keeps the evidence. The option stays opt-in.
- **Done (diagnostic)**: the step measured something and had no pass rule.
- **In progress**: the step runs now.
- **Next**: the step can start when the step before it is complete.
- **Blocked**: the step needs a decision or data. The row gives the reason.
- **Stopped**: the user stopped the line of work.

Counts are N at entrapment FDP_c 1 % on raw 01 / raw 02, unless the row says otherwise. Counts from q-values are decoy-count estimates.

## Now

| Step | WP | State | Record |
|---|---|---|---|
| Voxel pilot, 8 cells: no cell passes. Entrapment AUC 0.47–0.50 (base 0.516), no N gain over the voxel control | R-WP3 | Done (fail) | [RESULT_PILOT.md](gates/voxel-scoring-model/RESULT_PILOT.md) |
| Next direction for the voxel model | R-WP3 | Done: user chose heuristic features + semi-supervised rounds (10 November 2026) | [RESULT_PILOT.md](gates/voxel-scoring-model/RESULT_PILOT.md) |
| Voxel v2: 4 arms on raws 01/02: no arm passes. Best arm FS 921 / 1040 at 1 % (base 1024 / 1288), entrapment AUC above base | R-WP3 | Done (fail) | [RESULT_V2.md](gates/voxel-scoring-model/v2/RESULT_V2.md) |
| Voxel v3 step 0: brightness check not confirmed (drop ≤ 0.0011), voxel score as LDA feature does not help (V0 −; FS ±0) | R-WP3 | Done (fail) | [RESULT_STEP0.md](gates/voxel-scoring-model/v3-step0/RESULT_STEP0.md) |
| Voxel v3 phase 1: heuristic-only controls. LDA 1267 / 1261, MLP+BCE 1094 / 1100, linear+BCE 1300 / 1146. Any BCE model prefers entrapment (AUC 0.554–0.558 against 0.531): the objective is the bottleneck | R-WP3 | Done | [RESULT_P1.md](gates/voxel-scoring-model/v3-p1/RESULT_P1.md) |
| Voxel v3 phase 2a: objectives for one network (BCE all negatives, Fisher loss): none acceptable. Cause: positive-class covariance; decoy-only-covariance LDA keeps N (1259 / 1244) and AUC 0.527 | R-WP3 | Done (fail, cause found) | [RESULT_P2a.md](gates/voxel-scoring-model/v3-p1/RESULT_P2a.md) |
| Voxel v3 phase 2b gate 1: one network, decoy-whitened objective, 11 features (D11): PASS. 1348 / 1334 at 1 % (base 1024 / 1288), entrapment AUC 0.511 / 0.514 | R-WP3 | Done (pass) | [RESULT_P2b.md](gates/voxel-scoring-model/v3-p2b/RESULT_P2b.md) |
| Voxel v3 phase 2b gate 2: same objective with the voxel encoder (DV) | R-WP3 | **In progress** | [PLAN_V3_P2b.md](gates/voxel-scoring-model/v3-p2b/PLAN_V3_P2b.md) |
| C++ `OpenDIAScoringModel.{h,cpp}` and its own tool | R-WP3 | Blocked: needs a pilot cell that passes | [DESIGN.md](gates/voxel-scoring-model/DESIGN.md) |
| Full gate of `OpenDIAScoringModel` on raws 01/02 against the base (1029 / 1298) | R-WP3, R-WP5 | Blocked: needs the C++ model | |
| Pass rule: change from N@1 % to a paired bootstrap rule | R-WP4 | Blocked: needs user approval | [metric_noise.txt](gates/astral-native-byion/single-scorer/metric_noise.txt) |
| Raw 03 final check | R-WP6 | Blocked: needs user approval | |

## R-WP1: Streaming search

| Step | State | Record |
|---|---|---|
| Streaming sparse-grid search, multi-apex candidates, `-features_out` | Done | [R1](work-packages/R1-streaming-search.md) |
| Search inside calibrated windows (`-calibration`) | Done | [R1](work-packages/R1-streaming-search.md) |
| Native Thermo `.raw` input | Done | [R2 raw 01 section](work-packages/R2-calibration.md) |
| Raw-point export `-out_xipm` (OpenSWATH `.xipm` Parquet) | Done | [R4 post-mortem](work-packages/R4-error-control.md) |
| Consensus peak-group picker and anchor boundary labels (`e8cb5892b1`) | Done (diagnostic) | [peak-boundaries](gates/astral-raw01-peak-boundaries/RESULT.md) |
| Candidate recall at calibrated windows, by cause | Open | [R1](work-packages/R1-streaming-search.md) |

## R-WP2: Calibration

| Step | State | Record |
|---|---|---|
| `OpenDIACalibrate` on HeLa (1278 anchors) | Done (pass) | [R2](work-packages/R2-calibration.md) |
| Frozen design on Astral raw 01: 0 anchors | Done (fail) | [R2](work-packages/R2-calibration.md) |
| Windowed second pass `-window_panel_targets` (raw 01: 512 anchors) | Done (pass) | [R2](work-packages/R2-calibration.md) |
| Low-input runs (500 pg) | Open | [R2](work-packages/R2-calibration.md) |

## R-WP3: Scorer

| Step | State | Record |
|---|---|---|
| Native CNN crossfit, `-pooled_tdc` | Done | [R3](work-packages/R3-trace-cnn.md) |
| Stage-2 LDA (`-stage2`), `search` group (raw 01: 717 at q ≤ 0.01) | Done (pass) | [stage2-search](gates/astral-raw01-stage2-search/RESULT.md) |
| Stage-2 `isotope` group | Done (fail) | [stage2-isotope](gates/astral-raw01-stage2-isotope/RESULT.md) |
| Stage-2 `boundary` group (+4.0 %, rule was +5 %) | Done (fail) | [stage2-boundary](gates/astral-raw01-peak-boundaries/stage2-boundary/RESULT.md) |
| Stage-2 `byion` group: 1029 / 1298 (+36 % / +26 %). **Frozen native base.** | Done (pass on N) | [byion raw01](gates/astral-native-byion/raw01/RESULT.md), [raw02](gates/astral-native-byion/raw02/RESULT.md) |
| Stage-2 `shape` group: N up, q ≤ 0.05 calibration fails | Done (fail) | [shape01](gates/astral-native-byion/shape01/RESULT.md) |
| TOPAZ-style CNN arms A–D on raw 01 | Done (fail) | [cnn-topaz](gates/astral-raw01-cnn-topaz/RESULT.md) |
| Single TOPAZ scorer, rounds 1–5 (seeds, augmentation, decoy negatives, linear head, target negatives) | Done (fail) | [single-scorer](gates/astral-native-byion/single-scorer/RESULT.md), CALIB §13 |
| Round 6: transition refinement as features (R1 1028 / 1160) and as a mask (R2 520 / 679) | Done (fail) | [RESULT6.md](gates/astral-native-byion/single-scorer/RESULT6.md), CALIB §14 |
| Round 7: refined peak picking (P1 1024 / 1323, in the noise) | Done (fail) | [RESULT7.md](gates/astral-native-byion/single-scorer/RESULT7.md), CALIB §15 |
| Further work on CNN + LDA | **Stopped** (user decision, 10 November 2026) | |
| Single voxel model: design, training set, pilot rule (frozen before training) | Done | [DESIGN.md](gates/voxel-scoring-model/DESIGN.md) |
| Voxel pilot: export (`-out_xipm`, candidate table byte-identical to base) and voxel build | Done | [voxel-scoring-model/](gates/voxel-scoring-model/) |
| Voxel pilot: 8 cells (no cell passes, see RESULT_PILOT.md) | Done (fail) | [RESULT_PILOT.md](gates/voxel-scoring-model/RESULT_PILOT.md) |
| Next direction for the voxel model | Done: user chose heuristic features + semi-supervised rounds | [RESULT_PILOT.md](gates/voxel-scoring-model/RESULT_PILOT.md) |
| Voxel v2: 4 arms on raws 01/02: no arm passes. Best arm FS 921 / 1040 at 1 % (base 1024 / 1288), entrapment AUC above base | Done (fail) | [RESULT_V2.md](gates/voxel-scoring-model/v2/RESULT_V2.md) |
| Voxel v3 step 0: brightness check not confirmed, voxel score as LDA feature does not help | Done (fail) | [RESULT_STEP0.md](gates/voxel-scoring-model/v3-step0/RESULT_STEP0.md) |
| Voxel v3 phase 1: heuristic-only controls; the BCE objective is the bottleneck | Done | [RESULT_P1.md](gates/voxel-scoring-model/v3-p1/RESULT_P1.md) |
| Voxel v3 phase 2a: objectives for one network: none acceptable, cause found (positive-class covariance) | Done (fail, cause found) | [RESULT_P2a.md](gates/voxel-scoring-model/v3-p1/RESULT_P2a.md) |
| Voxel v3 phase 2b gate 1: one network, decoy-whitened objective, 11 features: PASS | Done (pass) | [RESULT_P2b.md](gates/voxel-scoring-model/v3-p2b/RESULT_P2b.md) |
| Voxel v3 phase 2b gate 2: same objective with the voxel encoder | **In progress** | [PLAN_V3_P2b.md](gates/voxel-scoring-model/v3-p2b/PLAN_V3_P2b.md) |
| C++ `OpenDIAScoringModel` | Blocked: needs pilot | |
| Full gate of `OpenDIAScoringModel` on raws 01/02 against the base (1029 / 1298) | Blocked: needs the C++ model | |

## R-WP4: Error control

| Step | State | Record |
|---|---|---|
| Family cross-fitting, paired TDC, pooled TDC | Done | [R4](work-packages/R4-error-control.md) |
| First entrapment result, raw 02, frozen design (FDP_c 0.028) | Done (fail) | [R4](work-packages/R4-error-control.md) |
| Decoy-q calibration diagnosis (decoy families, chimera, entrapment truth, commonness) | Done (diagnostic) | CALIB §1–8 |
| Label-free bulk-null density-ratio TDC (offline, both raws) | Done (pass) | CALIB §9 |
| Native bulk-null: needs a full HGB, not ported | Done (diagnostic) | CALIB §10 |
| Pipeline carriers of the 5 % excess: none found | Done (diagnostic) | CALIB §11 |
| Commonness-matched decoys, native | Done (fail) | CALIB §12 |
| Metric noise: bootstrap 90 % of N@1 % is about ±200 | Done (diagnostic) | [metric_noise.txt](gates/astral-native-byion/single-scorer/metric_noise.txt) |
| Calibrated q at q ≤ 0.05 (now 0.065 / 0.074) | Open | |
| Pass rule change to a paired bootstrap rule | Blocked: needs user approval | |

## R-WP5: ID gap

| Step | State | Record |
|---|---|---|
| HeLa loss funnel | Done (diagnostic) | [R5](work-packages/R5-id-gap.md) |
| Astral gap against DIA-NN and same-library OSW + Percolator | Done (diagnostic) | [diann-reference-gap](gates/astral-diann-reference-gap/GAP.md) |
| Raw 01 rethink: ceiling, oracle, RT lever | Done (diagnostic) | [REDESIGN.md](gates/astral-raw01-rethink/REDESIGN.md) |
| Raw 01/02 gap to OSW is score formation on the right peak | Done (diagnostic) | [gap2](gates/astral-raw01-rethink/gap2/) |
| Native base against OSW + Percolator: raw 02 1298 / 1287, raw 01 1029 / 1151 | Done | [byion raw01](gates/astral-native-byion/raw01/RESULT.md) |
| Close the raw 01 gap with the single voxel model | Blocked: needs pilot and C++ model | |

## R-WP6: Acceptance

| Step | State | Record |
|---|---|---|
| Raw 03 check under a new hashed gate | Blocked: needs user approval | |
| Native Bruker `.d` and Thermo `.raw` | Partial: both stream | [R6](work-packages/R6-acceptance.md) |
| Packaging (LibTorch optional, CTD/INI docs) | Open | [R6](work-packages/R6-acceptance.md) |
| Split the branch into reviewable PRs | Open | [R6](work-packages/R6-acceptance.md) |
