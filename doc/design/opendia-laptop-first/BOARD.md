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

Rows of 12-18 November 2026. Earlier rows (voxel v2-v19, CNN + LDA) are in the R-WP sections and in git history.

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
| Close the gap to DIA-NN. Astral (from .raw): 81-84 % at FDP_c 1 %. HeLa diaPASEF after RTCENTRE + IMGRID: 56 % (design) / 63 % (held out), 84 % of OpenSWATH SVM | R-WP5 | In progress | [RESULT_GAP.md](gates/extraction-parity/gap/RESULT_GAP.md), [RESULT_IMGRID.md](gates/single-binary/imgrid/RESULT_IMGRID.md) |
| HeLa diaPASEF final benchmark: OpenDIA 11,653 vs DIA-NN 28,950 (40 %), OpenSWATH SVM 19,132 (61 %); below the earlier rewrite (13,253), so the Astral changes regress diaPASEF | R-WP5, R-WP6 | Done (diagnostic) | [RESULT_HELA.md](gates/single-binary/bench/RESULT_HELA.md) |
| HeLa diagnosis: uncalibrated 1/K0 in the evidence pass (fixed 0e9e4fc3a5), early-gradient RT map bias (fixed by RTCENTRE 6dd893b954), mobility-blind grid/interference (fixed by IMGRID 963543c7f0); still open: b/y-ion floor, no 12-fragment HeLa library | R-WP2, R-WP5 | Done (diagnostic) | [DIAG_HELA.md](gates/single-binary/bench/DIAG_HELA.md) |
| Fair same-laptop benchmark (raw 01, 8 threads, both read .raw): OpenDIA 10:14, 15.6 GB, 172 % CPU; DIA-NN 7:23, 6.8 GB, 765 % CPU | R-WP6 | Done (diagnostic) | [README snapshot](README.md), [bench/](gates/single-binary/bench/) |
| Raw-pass speed: m/z-bucketed evidence pass (271 → 111 s) and level-aware Thermo read-ahead; laptop raw 01 10:14 → 7:18 (DIA-NN 7:23), byte-identical | R-WP6 | Done | [README snapshot](README.md) |
| Calibration memory (14.4 GB peak after eaa2c95cd0, from the 600k-member windowed pass; DIA-NN 6.8 GB) | R-WP6 | Open | [README snapshot](README.md) |
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
| HeLa regression diagnosis: 1/K0 bug, RT map and mobility-blind grid fixed (11,653 → 16,116); open: b/y-ion floor, 12-fragment HeLa library | Done (diagnostic) | [DIAG_HELA.md](gates/single-binary/bench/DIAG_HELA.md) |
| HeLa remaining gap (18 Nov): ~80 % right peak below cut; labelled ceiling on our 24 features 14,907 vs OpenSWATH sub-scores 17,822 (both 19,042); missing evidence (disjoint groups, corrected): xcorr co-elution +963, xcorr shape +739, fragment isotope +726, mass deviation +617, intensity/S/N +522; all together +1,816 | Done (diagnostic) | [GAP_HELA_20261118.md](gates/single-binary/bench/GAP_HELA_20261118.md) |
| OpenSWATH-style peak-group scores (PLAN_PGSCORES, 10 features): D 16,116 → 18,792 (+16.6 %), held-out H 24,000 → 28,020 (+16.8 %); Astral N@FDP_c 1 % +11.8 / +2.4 / +6.2 % (raw 03 held out); fails agreement on D (0.9104 < 0.9249) and FDP_c on raw 02 (0.0103 median, 0.0129 max); removed | Done (fail) | [RESULT_PGSCORES.md](gates/single-binary/pgscores/RESULT_PGSCORES.md) |
| Which peak-group score breaks q (PRECHECK_PGQ): none; every pg score's entrapment-vs-decoy AUC is 0.45-0.55; the raw 02 excess is real (pooled FDP_c 0.0100 vs CUR 0.0045) but absent on raws 01/03 (pooled over all three: 0.0072); no subset to drop | Done (diagnostic) | [PRECHECK_PGQ.md](gates/single-binary/pgscores/PRECHECK_PGQ.md) |
| Are the extra non-DIA-NN IDs real? (PLAUS_D): of NEW's 3,053 additions on D, 2,359 confirmed by DIA-NN (real); of the other 694, ~270-500 false (DIA-NN PEP and target-decoy bounds, Q6), so q is optimistic in the diaPASEF tail; `-pg_scores` restored as opt-in (a2aa8fc56c, default unchanged) | Done (diagnostic) | [PLAUS_D.md](gates/single-binary/pgscores/PLAUS_D.md) |
| Peak vs identity per score (Q6): no pg score is peak-only (decoys on a real peak AUC 0.48-0.53 against decoys elsewhere); identity AUCs 0.85-0.90 | Done (diagnostic) | [Q6_PEAK_IDENTITY.md](gates/single-binary/pgscores/Q6_PEAK_IDENTITY.md) |
| Power of the calibration criterion (Q7): raw 02 NEW 10 entrapments is what an exact 1 % run gives half the time (P = 0.51); per-raw FDP_c rule cannot tell calibrated from not; pooled raws 01-03 NEW 0.0073 vs CUR 0.0088; PRECHECK_PGQ corrected | Done (diagnostic) | [Q7_NEXT.md](gates/single-binary/pgscores/Q7_NEXT.md) |
| HeLa entrapment library (predicted, 58,998 peptide pairs, codon) + PLAN_PGENT: `-pg_scores` N@FDP_c 1 % D 9,701 -> 12,033 (+24 %), held-out H 23,136 -> 26,924 (+16 %); fails absolute FDP_c <= 0.012 (NEW 0.0166 / 0.0186, CUR 0.0176 / 0.0167); stays opt-in | Done (fail) | [RESULT_PGENT.md](gates/single-binary/pgent/RESULT_PGENT.md) |
| Shuffle-decoy precheck (D seed 0): calibration not fixed (CUR 1.63 -> 1.32 %, NEW 1.61 -> 1.85 %, rho 1.1-1.9 either way); excess is pipeline-wide, decoy-method swap is not the fix | Done (diagnostic) | [PRECHECK_SHUF.md](gates/single-binary/pgent/PRECHECK_SHUF.md) |
| (b) chosen. Foreign-proteome entrapment (Arabidopsis, D seed 0): the ~1.7 % excess is mostly near-homolog entrapments (accepted ones 3-4x enriched for Hamming <= 2 to human); post-hoc homolog-free FDP_c CUR 1.16 %, NEW 0.97 % (withdrawn by PGENT2); `-pg_scores` +22 % at 1 % true FDR | Done (diagnostic) | [PRECHECK_FOREIGN.md](gates/single-binary/pgent/PRECHECK_FOREIGN.md) |
| PLAN_PGENT2 (homolog-free Arabidopsis entrapment at construction): FDP_c at q <= 0.01 CUR 1.43 %, NEW 1.59 % (fail crit 1); N@FDP_c 1 % +17 %; PRECHECK_FOREIGN's 'near nominal' withdrawn; excess is real (~1.5x) | Done (fail) | [VERDICT_DESIGN.md](gates/single-binary/pgent2/VERDICT_DESIGN.md) |
| Entrapment bounds (BOUNDS): the reported FDP_c is the combined upper bound; lower bound N_E/N is 0.87-0.99 % in every diaPASEF run, so q is not shown to fail (true FDP ~0.9-1.9 %, inconclusive); r ~1.1-1.6 too small | Done (diagnostic) | [BOUNDS.md](gates/single-binary/pgent2/BOUNDS.md) |
| r ~7 homolog-free foreign entrapment (BOUNDS_R5): tight bounds show diaPASEF q <= 0.01 fails, true FDP 1.6-1.8 % (lower bound 1.58 %, 90 % CI from 1.37 %), same with and without `-pg_scores`; `-pg_scores` +51 % at 1 % true FDP; runs take ~1 h 50 min, 90 GB | Done (diagnostic) | [BOUNDS_R5.md](gates/single-binary/pgent2/BOUNDS_R5.md) |
| Locate the q asymmetry (b3): paired TDC fair overall but false targets beat their own pseudo-reverse decoy 63-64 % among accepted pairs (PAIRWIN); with decoys predicted like targets this drops to 50 % and paired-only q <= 0.01 is calibrated (lower bound 0.42-0.50 %, upper 0.84-0.99 %). Cause: decoys carry the target's predicted RT/IM/intensities (PRECHECK_PDEC). Also found: OpenMS PQP reader drops peptides without a gene mapping | Done (diagnostic) | [PRECHECK_PDEC.md](gates/single-binary/pgent2/PRECHECK_PDEC.md) |
| PLAN_PDEC (frozen 9feed18a6d), D verdict FAIL crit 1/2: DIALibGen-predicted decoys at r~7, q<=0.01 true FDP D 1.09-1.27 %, H 1.09-1.27 % (was 1.58-1.80 % pseudo-reverse); N@true 1 % D 8,318, H 20,932 (CUR) | Done (fail D; H reported) | [RESULT_H.md](gates/single-binary/pdec/RESULT_H.md) |
| OpenDIA-side decoy prediction: OpenDIAWorkflow -decoys predict (01d5ded04a, opt-in). Precheck 2 end to end (D s0, r 1.56): 99.9 % paired, q<=0.01 true FDP 0.62-1.02 % CUR, 0.53-0.86 % NEW, N@true 1 % about -3 % vs pseudo-reverse; 12-14 min, 29 GB. Frozen r~7 plan after the PLAN_PDEC verdict | Open (precheck met) | [DESIGN_OPENDIA_DECOYS.md](gates/single-binary/pdec/DESIGN_OPENDIA_DECOYS.md) |
| PLAN_ODP (frozen c3e71db6a9), D verdict FAIL crit 1/2 narrowly: OpenDIA -decoys predict at r~7, q<=0.01 true FDP D 1.06-1.21 %, H 1.06-1.26 %; same as the DIALibGen route; pair win on H 0.50 but the excess is unchanged, so the residual is not decoy asymmetry | Done (fail D; H reported) | [RESULT_H.md](gates/single-binary/pdec/RESULT_H.md) |
| PRECHECK_TDC (diagnostic): residual is not the q estimator; at q<=0.01 accepted entrapments outnumber their own accepted decoys 1.31-1.38x (1.7-1.9x at q<=0.005); accepted entrapments share >=3 fragments with a co-eluting identified HeLa target 2.7x as often as decoys (H 14.0 % vs 5.2 %). Next: interference-feature precheck (one -features_out run) | Done (diagnostic); next open | [PRECHECK_TDC.md](gates/single-binary/pdec/PRECHECK_TDC.md) |
| PRECHECK_INTERF (diagnostic): interference features flag the shared-signal false targets (53 % fully explained; AUC 0.74-0.75 vs decoys) but the scorer does not use them. A fair label-free rule (drop winners claimed by a stronger co-eluting winner, q recomputed) cuts the true FDP only ~0.05 pp at -0.3 to -0.9 % IDs (one seed); not enough alone | Done (diagnostic); decision open | [PRECHECK_INTERF.md](gates/single-binary/pdec/PRECHECK_INTERF.md) |
| PRECHECK_QCORR (diagnostic): a label-free shared-signal excess correction added to the TDC numerator did not move the q ≤ 0.01 cut on D or H seed 0, even at 20× strength; a shared-signal score penalty changed at most 0.2 % of IDs and did not improve the 1.11 / 1.27 % D or 1.07 / 1.22 % H bounds | Done (diagnostic); rejected | [PRECHECK_QCORR.md](gates/single-binary/pdec/PRECHECK_QCORR.md) |
| PRECHECK_IDGAP (diagnostic, r~7 D s0): label-free GBM 8,052 and a DIA-NN-labelled GBM ceiling 7,798 against the tool's 7,968 at FDP_c 1 %, so the gap is the features, not the scorer; the score is at chance in the lowest abundance quartile (AUC 0.514); 4 of 24 features are dead; a larger calibration panel restores anchors (140 to 1,521) but not the RT tolerance (41 to 38 s) | Done (diagnostic) | [PRECHECK_IDGAP.md](gates/single-binary/pdec/PRECHECK_IDGAP.md) |
| PLAN_EXTRA (frozen 5df717133a), `-extra_scores` (12 MS1 / mobility / mass-error features, 83a0f6c14a): D PASS all criteria, N@FDP_c 1 % 8,233 to 11,332 (+37.6 %, 44 % of DIA-NN), lower bound 1.06 to 1.16 %, combined 1.21 to 1.33 % (within tolerance), pair win 0.515. H (reported): 20,745 to 25,643 (+23.6 %), calibration unchanged (1.07 / 1.23 %), pair win 0.483. Now the default (969b474a51, byte-checked); `-no_extra_scores` restores the old scorer | Done (pass D; H reported; default) | [RESULT_EXTRA_D.md](gates/single-binary/pdec/RESULT_EXTRA_D.md), [RESULT_EXTRA_H.md](gates/single-binary/pdec/RESULT_EXTRA_H.md) |
| PRECHECK_IMPICK (diagnostic): mobility does not separate right from wrong candidates (RT does, 0.82); 16 % of DIA-NN IDs have no right candidate although inside both windows, so 2-D picking targets candidate recall | Done (diagnostic) | [PRECHECK_IMPICK.md](gates/single-binary/pdec/PRECHECK_IMPICK.md) |
| RT x IM candidate picking `-im_peaks` (07131138b8, opt-in): r~1.56 s0 right candidates 84 to 94 %, +3.5 %; PLAN_IM2D (frozen 38273668e7) D FAIL crit 1: N@FDP_c 1 % 11,332 to 11,420 (+0.8 %, needs +3 %), calibration and pair win unchanged, +49 % wall, +13 GB. Recovered peaks still score low. H (reported): 25,643 to 26,302 (+2.6 %), calibration unchanged, +47 % wall. Stays opt-in | Done (fail D; H reported) | [RESULT_IM2D_D.md](gates/single-binary/pdec/RESULT_IM2D_D.md), [RESULT_IM2D_H.md](gates/single-binary/pdec/RESULT_IM2D_H.md), [PRECHECK_IM2D.md](gates/single-binary/pdec/PRECHECK_IM2D.md) |
| PRECHECK_SCORER36 (diagnostic, r~7 D s0, 36-feature default): tool 11,258 at FDP_c 1 %, label-free GBM 11,161, 2x GBM 11,368, DIA-NN-labelled ceiling 11,122, so the scorer is still not the limit; of 25,881 DIA-NN IDs, 4,605 (17.8 %) have no window within 6 s and ~10,000 score below the decoy tail (Q1 TPR 0.22 at decoy FPR 1e-4). Side lead: hard-negative label-free training +5 % (seed 0, calibration unchanged), needs its own precheck + plan | Done (diagnostic) | [PRECHECK_SCORER36.md](gates/single-binary/pdec/PRECHECK_SCORER36.md) |
| PRECHECK_NOCAND (diagnostic): 88 % of the 4,605 r~7 DIA-NN IDs with no candidate are merged into a neighbouring peak by the search valley rule (right peak present, 6/6 fragments, ~0.6 x the neighbour). Advanced `-valley_fraction`/`-max_peaks` (822ef8e362, defaults byte-identical): valley 0.8 + max 8 at r~1.56 s0 gives a right candidate for 84 to 98 % of DIA-NN IDs but only +1.8 % N@FDP_c 1 % at 1.9x wall, so not planned. Hard-negative training: +1.2 % then worse calibration at r~1.56, dropped. Next lever: faint-peak evidence | Done (diagnostic) | [PRECHECK_NOCAND.md](gates/single-binary/pdec/PRECHECK_NOCAND.md) |
| PRECHECK_FAINT (diagnostic, r~1.56 D s0, valley 0.8/max 8): 7 per-scan trace features (library cosine at apex, local S/N, ...) raise a DIA-NN-labelled ceiling TPR at decoy FPR 1e-3 0.509 to 0.535 and at 1e-2 0.658 to 0.693 (bootstrap 90 % excludes 0, both GBM seeds), i.e. +5-7 % relative, mostly ps_cos_apex and local_sn; Q1 still 0.22. Cheap to add, not a parity route | Done (diagnostic) | [PRECHECK_FAINT.md](gates/single-binary/pdec/PRECHECK_FAINT.md) |
| Narrower 1/K0 window `-im_tolerance_scale` 0.6, now default (OpenMS 9ce3370ee3, byte-identical to tested arm): PLAN_IMNARROW (frozen 07a4245472) D PASS (ab765855cc) N@FDP_c 1 % 11,332 to 11,815 (+4.3 %); H (reported) 25,643 to 26,647 (+3.9 %); every NEW seed above every CUR seed on both; calibration unchanged; no extra cost. 45 % / 65 % of DIA-NN | Done | [RESULT_IMNARROW_D.md](gates/single-binary/pdec/RESULT_IMNARROW_D.md), [RESULT_IMNARROW_H.md](gates/single-binary/pdec/RESULT_IMNARROW_H.md), [PRECHECK_IMNARROW.md](gates/single-binary/pdec/PRECHECK_IMNARROW.md), [PLAN_IMNARROW.md](gates/single-binary/pdec/PLAN_IMNARROW.md) |
| PRECHECK_NNCEIL (diagnostic, faint sample): a label-free 1-D trace CNN as one extra input adds as much as the 7 faint features (ceiling TPR @1e-2 0.66 to 0.70, Q1 +0.04-0.05) and overlaps them almost fully; alone it is far weaker than the hand features. Next net test needs mobility-resolved input | Done (diagnostic) | [PRECHECK_NNCEIL.md](gates/single-binary/pdec/PRECHECK_NNCEIL.md) |
| PRECHECK_IMFEAT (diagnostic, r~1.56 s0 on the narrow-window default): own-mobility-apex features (band6) add +2.3/+2.4 TPR pts at decoy FPR 3e-3 over the 36 (ceiling), Q1 +4 pts: proceed. faint2 nothing new at the narrow window. 2-D trace CNN alone 0.54 vs 1-D 0.26; stacked +1.5-1.8 (below bar). C++-feasible core-window 1/K0 histogram (core5) keeps the gain: +2.2/+2.2, proceed; tensor-only variant (tensor5) +1.8/+2.0, not taken. core5 implemented as opt-in `-mobility_apex` (OpenMS 8f113de11a, default byte-identical). Workflow check r~1.56 D 2 seeds: N@FDP_c 1 % +4 / -51 (inside seed spread), pooled combined bound 1.15 to 1.40 %, no decoy leak (AUC 0.49-0.50): ceiling gain not realised, no plan; flag kept opt-in | Done (diagnostic) | [PRECHECK_IMFEAT.md](gates/single-binary/pdec/PRECHECK_IMFEAT.md) |
| PQP reader dropped peptides without a gene mapping when GENE has rows; fixed (row order unchanged), class test (2982bd647d) | Done | [PRECHECK_PDEC.md](gates/single-binary/pgent2/PRECHECK_PDEC.md) |
| GAPDIAG (diagnostic, r~7 D s0 vs DIA-NN): of DIA-NN's 25,881 IDs we find 45.7 %; 34.0 % have our best window on the right peak but score low (co-elution near decoy level, our area vs DIA-NN quantity Spearman 0.54 vs 0.78); 8.2 % wrong peak; 12.2 % no window. RT and IM not the cause. Lead: 0.02 Da trace window (22-58 ppm) vs 8.3 ppm calibrated. Opt-in `-trace_tolerance_ppm` (OpenMS e19261618f, default byte-identical): r~1.56 D 10 ppm N@FDP_c 1 % +11.7 % / +9.8 % (2 seeds), combined bound 1.15 to 1.56 % | Done (diagnostic) | [GAPDIAG.md](gates/single-binary/pdec/GAPDIAG.md) |
| PLAN_TRPPM (frozen 0d80832844, sha 8c909603): `-trace_tolerance_ppm 10`, r~7 D seeds 0-4. **D PASS**: N@FDP_c 1 % 11,815 to 13,985 (+18.4 %, every NEW seed above every CUR seed), lower 1.19 to 1.10 %, combined 1.36 to 1.26 %, pair win 0.499, wall 0.95x, 34 GiB. 45 % to 54 % of DIA-NN on D. H (reported after the D verdict): 26,647 to 29,286 (+9.9 %, every NEW seed above every CUR seed), lower 1.06 to 1.15 %, combined 1.21 to 1.31 %; 65 % to 71 % of DIA-NN. Default 10 ppm in OpenMS ffc7daaa99 (byte-identical to tested arm; 0 = previous default byte-identical; 25/25) | Done (pass) | [RESULT_TRPPM_D.md](gates/single-binary/pdec/RESULT_TRPPM_D.md), [RESULT_TRPPM_H.md](gates/single-binary/pdec/RESULT_TRPPM_H.md), [PLAN_TRPPM.md](gates/single-binary/pdec/PLAN_TRPPM.md) |
| Memory rounds 1-4 (output-identical, 25/25): chunked grid pool, release library copy, Sample ids out, streaming TraceEvidence (OpenMS 85c794b211, d0bdf24540, f3f1fc046f, 0a49d8aa8e). r~1.56 D 11.28 to about 10.2-10.5 GB; r~7 D round 4 35.57 to 34.11 GB. Peak is extraction: tensors 12.2 GB, library 7 GB, grids 6.0 GB, IR windows 3.9 GB. No exact win left; DIA-NN's <7 GB needs an output-changing or structural change (own plan). Round 4 same-node A/B: wall +1.3 %, peak -1.5 GB | Done | [MEM_ROUNDS.md](gates/single-binary/workflow/MEM_ROUNDS.md) |
| PLAN_FAINT (frozen 20 Nov 2026, source `11117425cc`): opt-in `-faint_features` adds pooled trace cosine-at-apex and local-S/N features. r~1.56 precheck +6.3 % median N. r~7 D seeds 0-4: 13,985 to 14,549 (+4.0 %, below frozen +5 % criterion), lower 1.10 to 1.13 %, combined 1.26 to 1.29 %, pair win 0.499 to 0.512, cost passes. H not run; flag remains opt-in | Done (fail D; H not run) | [RESULT_FAINT_D.md](gates/single-binary/pdec/RESULT_FAINT_D.md), [PLAN_FAINT.md](gates/single-binary/pdec/PLAN_FAINT.md) |
| PRECHECK_MF (diagnostic, faint2 sample on the current default): fragment-level matched filter (per-fragment robust S/N, leave-one-out shape, fragment dominance) passes its pre-set gate. Raw-point mf: ceiling TPR @3e-3 +0.066/+0.071, label-free proxy on missed DIA-NN IDs 0.22 to 0.33-0.37; tensor-only mt: +0.034/+0.037 and +0.091/+0.103. C++ opt-in `-matched_filter` (OpenMS 34c36861f2, default byte-identical on r~1.56 D s0/s1, 25/25). r~1.56 D workflow check: N@FDP_c 1 % 14,908 to 16,243 (+9.0 %, both seeds), bounds 0.82/1.56 to 0.84/1.58 %, pair win 0.514, same wall | Done (pass) | [PRECHECK_MF.md](gates/single-binary/pdec/PRECHECK_MF.md) |
| PLAN_MT (frozen ec28ba922d): `-matched_filter` on r~7 D seeds 0-4. **D PASS**: N@FDP_c 1 % 13,985 to 14,895 (+6.5 %, every MT seed above every CUR seed), lower 1.10 to 1.20 %, combined 1.26 to 1.37 % (inside +0.2/+0.3 pp limits), pair win 0.533, wall 0.97x, 34.4 GB. 57 % of DIA-NN on D. H (reported after the D verdict): 29,286 to 29,873 (+2.0 %, every MT seed above every CUR seed), bounds unchanged 1.15/1.31 to 1.15/1.32 %, pair win 0.517; 72 % of DIA-NN. Now the default (OpenMS 01c5891ce0, byte-identical to the tested arm; `-no_matched_filter` byte-identical to the previous default; 25/25) | Done (pass) | [RESULT_MT_D.md](gates/single-binary/pdec/RESULT_MT_D.md), [RESULT_MT_H.md](gates/single-binary/pdec/RESULT_MT_H.md), [PLAN_MT.md](gates/single-binary/pdec/PLAN_MT.md) |
| PRECHECK_MF2 (diagnostic, over the new default): raw-point info adds little once mt is in. Full raw mf +0.018/+0.034 tool-missed; mass weighting alone nothing; own-mobility centring +0.017/+0.015. Implementable variants (ppm-weighted tensor, 16-bin IM tensor) fail the pre-set gate. No point-cloud implementation; next levers are H's RT calibration (p95 residual 52 s vs 8 s on D) and wrong-peak / no-candidate classes | Done (fail) | [PRECHECK_MF2.md](gates/single-binary/pdec/PRECHECK_MF2.md) |
| PRECHECK_RTIM (diagnostic): H's RT map is calibrated on H itself and is as good as D's relative to gradient length (p95 10.0 vs 9.8 % of span; both near the monotone floor). Sage-style per-run ridge on sequence composition, trained on the run's own q<=0.01 IDs (cross-fitted): **RT passes**, p95 on missed DIA-NN IDs D 22.6 to 16.4 s (-27 %), H 102.2 to 63.1 s (-38 %), below the monotone floor; ~3,000 IDs reach the full gain. **1/K0 fails** (-5 to -8 %). Next: opt-in second-pass RT correction in C++ | Done (RT pass, IM fail) | [PRECHECK_RTIM.md](gates/single-binary/pdec/PRECHECK_RTIM.md) |
| `-rt_refine` (OpenMS f115760e89) and `-rt_refine_anchors` (43da5895b1), opt-in, default byte-identical, 26/26. r~1.56 D 2-seed check: `-rt_refine` N@FDP_c 1 % 16,243 to 16,858 (+3.8 %, both seeds above), bounds 0.84/1.58 to 0.87/1.64 %, RT window 35.2 to ~26.3 s, wall 1.7-1.9x; anchors-only +1.9 % with a higher lower bound (not taken). PLAN_RTR (frozen 21 Nov) on r~7 D seeds 0-4, +3 % bar, H after the D verdict. **D PASS** (narrow): N@FDP_c 1 % 14,895 to 15,369 (+3.2 %, bar 15,342; seed 2 flat), calibration improves (lower 1.20 to 1.12 %, combined 1.37 to 1.27 %), pair win 0.527, wall 1.83x, 34.7 GB; RT window 41.4 to ~31.6 s. 59 % of DIA-NN on D. H running | Done (pass D); H running | [RESULT_RTR_D.md](gates/single-binary/pdec/RESULT_RTR_D.md), [PLAN_RTR.md](gates/single-binary/pdec/PLAN_RTR.md) |
| Near-homolog (sibling) competition: run-derived, via the interference graph | Open (idea) | [PRECHECK_FOREIGN.md](gates/single-binary/pgent/PRECHECK_FOREIGN.md) |
| b/y-ion floor; 12-fragment HeLa library | Open | [GAP_HELA_20261118.md](gates/single-binary/bench/GAP_HELA_20261118.md) |
| One-pass stage 3 centred on the stage-1 map (PLAN_RTCENTRE): HeLa 13,127 → 14,875 (+13.3 %), held-out HeLa 19,725 → 21,682 (+9.9 %), agreement with DIA-NN held, Astral identical; kept (6dd893b954) | Done (pass) | [RESULT_RTCENTRE.md](gates/single-binary/rtcentre/RESULT_RTCENTRE.md) |
| Mobility-resolved 12-fragment grid and interference (PLAN_IMGRID): HeLa 14,875 → 16,116 (+8.3 %), held-out 21,682 → 24,000 (+10.7 %), agreement held, Astral identical; kept (963543c7f0) | Done (pass) | [RESULT_IMGRID.md](gates/single-binary/imgrid/RESULT_IMGRID.md) |
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
| Fair same-laptop benchmark (raw 01, 8 threads, both read .raw): OpenDIA 10:14, 15.6 GB, 172 % CPU; DIA-NN 7:23, 6.8 GB, 765 % CPU | Done | [bench/](gates/single-binary/bench/) |
| Raw-pass speed: m/z-bucketed evidence pass (271 → 111 s) and level-aware Thermo read-ahead; laptop raw 01 10:14 → 7:18 (DIA-NN 7:23), byte-identical | Done | [README snapshot](README.md) |
| Calibration memory (14.4 GB peak after eaa2c95cd0, from the 600k-member windowed pass; DIA-NN 6.8 GB) | Open | [README snapshot](README.md) |
| HeLa diaPASEF benchmark: 11,653 IDs, 40 % of DIA-NN, 61 % of OpenSWATH SVM, below the 2026-11-07 rewrite (13,253) | Done (diagnostic) | [RESULT_HELA.md](gates/single-binary/bench/RESULT_HELA.md) |
| Native Bruker `.d` and Thermo `.raw` | Partial: both stream | [R6](work-packages/R6-acceptance.md) |
| Packaging (LibTorch optional, CTD/INI docs) | Partial: LibTorch off the default path | [R6](work-packages/R6-acceptance.md) |
| Split the branch into reviewable PRs | Open | [R6](work-packages/R6-acceptance.md) |
| Thermo RAW level-aware read-ahead split out: PR OpenMS#10925, approved by Timo; review comments addressed in abc2eb24af (transform MS2 4 threads on Astral: 25.4 s develop → 16.7 s) | Done (PR open, approved) | [#10925](https://github.com/OpenMS/OpenMS/pull/10925) |
| One-command `OpenDIAWorkflow`, deterministic, fail-closed CLI tests (`9d6bf16f10`) | Done | [RESULT_WORKFLOW.md](gates/single-binary/workflow/RESULT_WORKFLOW.md) |
