# OpenDIA rewrite: step board

[Master plan](README.md)

**Status, 16 November 2026.** This file lists the steps of each work package and the state of each step. Each row links to the record of the step. The script [`board_sync.py`](board_sync.py) copies the R-WP sections of this file to the GitHub project [OpenMS/38](https://github.com/orgs/OpenMS/projects/38). Edit this file first, then run the script. The "Now" table repeats rows from the R-WP sections, and the script does not copy it.

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

Rows of 12-16 November 2026. Earlier rows (voxel v2-v19, CNN + LDA) are in the R-WP sections and in git history.

| Step | WP | State | Record |
|---|---|---|---|
| Native-scan XIC IDs (PLAN_XIC_IDS): N up +390..+743 / +546..+1051 at 1 %, FDP_c at q ≤ 0.01 0.016 / 0.017 | R-WP3, R-WP4 | Done (fail, rule 1c) | [RESULT_XIC_IDS.md](gates/extraction-parity/RESULT_XIC_IDS.md) |
| Fragment-deficit features (PLAN_INTERF): raw 01 FDP_c 0.023 | R-WP3 | Done (fail) | [RESULT_INTERF.md](gates/extraction-parity/RESULT_INTERF.md) |
| Parent-complete library (PLAN_PAIRLIB): excess stays, FDP_c 0.032 / 0.017 | R-WP4 | Done (fail) | [RESULT_PAIRLIB.md](gates/extraction-parity/pairlib/RESULT_PAIRLIB.md) |
| Predicted decoys (PLAN_DECOYGEN, `-predict_decoys`): rule 1c passes in every arm; C++ arm E reproduces arm D | R-WP4 | Done (pass, arm E) | [RESULT_DECOYGEN.md](gates/extraction-parity/decoygen/RESULT_DECOYGEN.md) |
| Gap funnel against DIA-NN: losses are right-peak-below-cut (626 / 644) | R-WP5 | Done (diagnostic) | [RESULT_GAP.md](gates/extraction-parity/gap/RESULT_GAP.md) |
| PLAN_GAP: shared-signal interference (I) passes, 1694 / 1628 → 1879 / 1921; A, G fail. C++ OpenDIAInterference reproduces it | R-WP3, R-WP5 | Done (pass, I) | [RESULT_GAP_ARMS.md](gates/extraction-parity/gap/arms/RESULT_GAP_ARMS.md) |
| Dead-code removal (TOPAZ, MIL, stage-2 isotope/shape/refine, OpenDIARerank/Traces); outputs byte-identical | R-WP6 | Done | `86d82dfdd4` |
| C++ CrossFitScorer (D11 network, self-training, W rule, tail, pooled TDC): 98.6 / 99.0 % of Python median | R-WP3 | Done (pass) | [RESULT_STEP1.md](gates/single-binary/step1/RESULT_STEP1.md) |
| OpenDIAWorkflow, one binary on .raw: 1819 / 1870 (reused cal.), raw 02 at 97.1 % of Python (rule 98 %) | R-WP3, R-WP6 | Done (raw 01 pass, raw 02 fail) | [RESULT_WORKFLOW.md](gates/single-binary/workflow/RESULT_WORKFLOW.md) |
| Validation: determinism, report, run and global q recompute, fail-closed CLI tests | R-WP6 | Done (pass) | [RESULT_WORKFLOW.md](gates/single-binary/workflow/RESULT_WORKFLOW.md) |
| Self-calibrated RT spread: CALSPREAD fails offline; CALREFINE / CALREFINE2 (5 seeds) +5.8 % / −3.1 %, fail | R-WP2 | Done (fail) | [calspread](gates/single-binary/calspread/RESULT_CALSPREAD.md), [calrefine](gates/single-binary/calrefine/RESULT_CALREFINE2.md) |
| RT window width rule: density crossing undefined, P_win against the evidence; cross-fitted width on IDs (PLAN_RTWIN2) −2.2 % raw 01, held-out raw 03 −2.3 %, removed | R-WP2 | Done (fail) | [RESULT_RTWIN2.md](gates/single-binary/rtwin/RESULT_RTWIN2.md) |
| ID-based recalibration with RT-dependent window (PLAN_RTPROFILE): ends recovered, middle loses; raw 02 −1.05 %, raw 01 s4 FDP_c 0.0104, raw 03 −0.6 %, removed | R-WP2 | Done (fail) | [RESULT_RTPROFILE.md](gates/single-binary/rtwin/RESULT_RTPROFILE.md) |
| Close the gap to DIA-NN (DIA-NN from .raw: 2204 / 2356 / 2259 at FDP_c 1 %; OpenDIA 81-84 %) | R-WP5 | In progress | [RESULT_GAP.md](gates/extraction-parity/gap/RESULT_GAP.md) |
| Fair same-laptop benchmark (raw 01, 8 threads, both read .raw): OpenDIA 10:14, 15.6 GB, 172 % CPU; DIA-NN 7:23, 6.8 GB, 765 % CPU | R-WP6 | Done (diagnostic) | [README snapshot](README.md), [bench/](gates/single-binary/bench/) |
| Parallelise the single-threaded raw passes (calibration x2, search, evidence pass) and cut calibration memory | R-WP6 | Open | [README snapshot](README.md) |
| Status snapshot vs DIA-NN: 87 % / 84 % of DIA-NN at FDP_c 1 %; 67-72 % of DIA-NN's q ≤ 0.01 IDs shared, 93 % of ours in DIA-NN; 6-7.5 min, 22 GB per raw | R-WP5, R-WP6 | Done (diagnostic) | [README snapshot](README.md), [STATUS_vs_DIANN.txt](gates/single-binary/STATUS_vs_DIANN.txt) |
| Joint (presence-weighted) interference fit: no room, block I's explaining-away comes from absent claimants; labelled ceiling on current inputs about +100 / raw | R-WP5 | Done (diagnostic, not pursued) | [PRECHECK_JOINT.md](gates/extraction-parity/gap/PRECHECK_JOINT.md) |
| Scorer capacity (option 1): interaction inputs and self-trained trees fail offline; seed ensemble (PLAN_ENSEMBLE, m = 3) +0.2 % / +0.6 % (rule +1 %), raw 03 +0.3 %, removed | R-WP3, R-WP5 | Done (fail) | [RESULT_ENSEMBLE.md](gates/single-binary/ensemble/RESULT_ENSEMBLE.md) |
| Raw 03 held-out check: self-calibrating, `library12_full`, 1907 median, FDP_c ≤ 0.0096 | R-WP6 | Done | [RESULT_RTWIN2.md](gates/single-binary/rtwin/RESULT_RTWIN2.md) |

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
| Selection-free RT spread fit (`20f9f12335`) and RT-likelihood stage-3 ranking (`0e4a234762`) | Done | [RESULT_CALIB.md](gates/extraction-parity/decoygen/calib/RESULT_CALIB.md) |
| Self-calibrated RT window on the predicted-decoy library: spread inflated (frame mismatch), fixes fail (CALSPREAD, CALREFINE, CALREFINE2) | Done (fail) | [RESULT_CALREFINE2.md](gates/single-binary/calrefine/RESULT_CALREFINE2.md) |
| RT window width rule (PRECHECK_CROSSING, PLAN_RTWIN2 cross-fitted selection) | Done (fail, removed) | [RESULT_RTWIN2.md](gates/single-binary/rtwin/RESULT_RTWIN2.md) |
| ID-based recalibration, RT-dependent window (PLAN_RTPROFILE) | Done (fail, removed) | [RESULT_RTPROFILE.md](gates/single-binary/rtwin/RESULT_RTPROFILE.md) |

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
| Voxel v3 phase 2b gate 2: voxel encoder + 11 features: FAIL (voxel adds no IDs over D11) | Done (fail) | [RESULT_P2b.md](gates/voxel-scoring-model/v3-p2b/RESULT_P2b.md) |
| Voxel v3-seq: ReDeeM-style sequence encoder (DS11, DVS): FAIL (no gain at 2 %, no sequence tell) | Done (fail) | [RESULT_SEQ.md](gates/voxel-scoring-model/v3-seq/RESULT_SEQ.md) |
| Next step after v3-seq | Blocked: needs user decision | [RESULT_SEQ.md](gates/voxel-scoring-model/v3-seq/RESULT_SEQ.md) |
| Voxel v4 phase 1: ranking term and MS1 features: FAIL (no gain over D11 at 2 %) | Done (fail) | [RESULT_V4P1.md](gates/voxel-scoring-model/v4/RESULT_V4P1.md) |
| Voxel v4 phase 2: 12-fragment features: D11F and D11FW PASS (gain at 2 %, 1 % unstable on raw 02) | Done (pass) | [RESULT_V4P2.md](gates/voxel-scoring-model/v4/RESULT_V4P2.md) |
| C++ `OpenDIAScoringModel` | Blocked: needs pilot | |
| Full gate of `OpenDIAScoringModel` on raws 01/02 against the base (1029 / 1298) | Blocked: needs the C++ model | |
| Voxel v13: precursor-aware full-spectrum NNLS diagnostic on raw 01: D-vs-decoy ablation AUC 0.687, but entrapment-vs-decoy AUC 0.535 misses frozen ≤0.53 rule. No training; keep D11FWT and stop this variant | Done (fail) | [RESULT_V13.md](gates/voxel-scoring-model/v13/RESULT_V13.md) |
| Voxel v14b: Stage-D-like local-flank attenuation diagnostic, raw01: median score barely changes, only 2.64% of high-confidence donors enter observed Stage-D score range. F12 inputs largely discard amplitude; no augmentation model trained | Done (diagnostic) | [RESULT_V14b.md](gates/voxel-scoring-model/v14/RESULT_V14b.md) |
| Shared-signal interference features, C++ `OpenDIAInterference` (`6560d83cb2`) | Done (pass) | [RESULT_GAP_ARMS.md](gates/extraction-parity/gap/arms/RESULT_GAP_ARMS.md) |
| C++ `CrossFitScorer` (`b851c0504a`) replaces the Python D11 trainer | Done (pass) | [RESULT_STEP1.md](gates/single-binary/step1/RESULT_STEP1.md) |
| `OpenDIAWorkflow` single binary (`794ee1565e`): raw 01 pass, raw 02 97.1 % (rule 98 %) | Done (partial) | [RESULT_WORKFLOW.md](gates/single-binary/workflow/RESULT_WORKFLOW.md) |

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
| Predicted decoys (`-predict_decoys`) remove the entrapment excess (rule 1c passes) | Done (pass) | [RESULT_DECOYGEN.md](gates/extraction-parity/decoygen/RESULT_DECOYGEN.md) |

## R-WP5: ID gap

| Step | State | Record |
|---|---|---|
| HeLa loss funnel | Done (diagnostic) | [R5](work-packages/R5-id-gap.md) |
| Astral gap against DIA-NN and same-library OSW + Percolator | Done (diagnostic) | [diann-reference-gap](gates/astral-diann-reference-gap/GAP.md) |
| Raw 01 rethink: ceiling, oracle, RT lever | Done (diagnostic) | [REDESIGN.md](gates/astral-raw01-rethink/REDESIGN.md) |
| Raw 01/02 gap to OSW is score formation on the right peak | Done (diagnostic) | [gap2](gates/astral-raw01-rethink/gap2/) |
| Native base against OSW + Percolator: raw 02 1298 / 1287, raw 01 1029 / 1151 | Done | [byion raw01](gates/astral-native-byion/raw01/RESULT.md) |
| Close the raw 01 gap with the single voxel model | Blocked: needs pilot and C++ model | |
| Gap funnel against DIA-NN (arm E): right peak below cut dominates | Done (diagnostic) | [RESULT_GAP.md](gates/extraction-parity/gap/RESULT_GAP.md) |
| Close the remaining 12-21 % gap to DIA-NN | In progress | [README snapshot](README.md) |
| Joint interference fit precheck: no room (absent claimants carry block I) | Done (diagnostic, not pursued) | [PRECHECK_JOINT.md](gates/extraction-parity/gap/PRECHECK_JOINT.md) |
| Scorer seed ensemble (PLAN_ENSEMBLE): +0.2 % / +0.6 %, removed | Done (fail) | [RESULT_ENSEMBLE.md](gates/single-binary/ensemble/RESULT_ENSEMBLE.md) |
| Wrong / missed-peak losses precheck: ~100 / raw, no cheap fix (RT window ends, weak peaks behind the 5-peak cap, W rule already better than max-score) | Done (diagnostic, not pursued) | [PRECHECK_PEAK.md](gates/extraction-parity/gap/PRECHECK_PEAK.md) |
| Run-specific fragment choice precheck: re-selection AUC 0.708→0.710 / 0.725→0.721, no room | Done (diagnostic, not pursued) | [PRECHECK_FRAG.md](gates/single-binary/frag/PRECHECK_FRAG.md) |
| MS1 isotope evidence (PLAN_MS1): N@FDP_c 1 % +5.7 % / +2.1 % on raws 01/02, fails q ≤ 0.01 and specificity rules on raw 01; held-out raw 03 −0.5 %; removed | Done (fail) | [RESULT_MS1.md](gates/single-binary/ms1/RESULT_MS1.md) |

## R-WP6: Acceptance

| Step | State | Record |
|---|---|---|
| Raw 03 held-out check under each new hashed plan | Done for PLAN_RTWIN2 / PLAN_RTPROFILE (current setup 1907 median) | [RESULT_RTWIN2.md](gates/single-binary/rtwin/RESULT_RTWIN2.md) |
| Memory: the 22 GB peak came from the windowed calibration pass, not the scoring library. Trimmed calibration traces: 21.6 → 15.6 GB, outputs byte-identical (a132e9104e) | Done | [README snapshot](README.md) |
| Native Bruker `.d` and Thermo `.raw` | Partial: both stream | [R6](work-packages/R6-acceptance.md) |
| Packaging (LibTorch optional, CTD/INI docs) | Partial: LibTorch off the default path | [R6](work-packages/R6-acceptance.md) |
| Split the branch into reviewable PRs | Open | [R6](work-packages/R6-acceptance.md) |
| One-command `OpenDIAWorkflow`, deterministic, fail-closed CLI tests (`9d6bf16f10`) | Done | [RESULT_WORKFLOW.md](gates/single-binary/workflow/RESULT_WORKFLOW.md) |
