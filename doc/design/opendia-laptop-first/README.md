# OpenDIA rewrite: master plan

**Status, 10 November 2026.** This plan replaced the laptop-first WP-00 to WP-10 plan on 7 November 2026. The old plan is closed and kept unchanged under [archive/legacy-wp-plan-20261107](archive/legacy-wp-plan-20261107/README.md).

The rewrite lives on [`OpenMS/OpenMS:feature/opendia-clean-rewrite`](https://github.com/OpenMS/OpenMS/tree/feature/opendia-clean-rewrite).
All of its code is under `src/topp/OpenDIA/`, and it makes no changes to the shared OpenMS API. The branch tip is `dfbe38b470`.
It is an experimental branch, not a pull request. It must not merge into `develop` until the acceptance gates of R-WP6 pass.
Detailed evidence lives in the branch's [`src/topp/OpenDIA/README_TRACE_MODEL.md`](https://github.com/OpenMS/OpenMS/blob/feature/opendia-clean-rewrite/src/topp/OpenDIA/README_TRACE_MODEL.md) and in [`gates/`](gates/).

## Current state, 10 November 2026

**Goal.** OpenDIA must identify at least as many precursors as OpenSWATH + Percolator, with honest error control and no teacher.

**Design data.** ProteoBench Astral entrapment raws 01 and 02 are the design runs. Raw 03 is held back for the final check. Raw 03 needs the user's approval.

**Metric.** N is the number of target precursors at an entrapment FDP_c of 1 %, with FDP_c = N_e (1 + 1/r) / N and r = 148983 / 151239. The bootstrap 90 % interval of N@1 % is about ±200 (`gates/astral-native-byion/single-scorer/metric_noise.txt`). A change of the pass rule to a paired bootstrap rule is proposed. It is not approved.

**Native base (frozen).** OpenDIACalibrate, then the OpenDIA search, then `OpenDIACNN -stage2` with the feature set `search+boundary+byion` (CNN logit plus a stage-2 LDA).

| N at entrapment FDP_c 1 % | raw 01 | raw 02 |
|---|---|---|
| Native base (CNN + LDA, `search+boundary+byion`) | 1029 | 1298 |
| OSW + Percolator, same library (reference only) | 1151 | 1287 |

- On raw 02, the native base is at the OSW + Percolator level. On raw 01, it is at 89 %.
- Plain pooled q at q ≤ 0.01 controls the error on both raws (FDP_c 0.013 and 0.007).
- At q ≤ 0.05, plain pooled q is anti-conservative (FDP_c 0.065 and 0.074). This is open in R-WP4.

**Direction (user decision, 10 November 2026).** Stop the work on CNN + LDA. Replace both with **one voxel scoring model** in new files, `OpenDIAScoringModel.{h,cpp}`, with its own tool. Train the model with decoy-like target negatives, so that it separates false targets from true targets. The design and the frozen pilot rule are in [`gates/voxel-scoring-model/`](gates/voxel-scoring-model/DESIGN.md).

## Why the direction changed

The legacy plan kept OpenSWATH's chromatogram and feature-scoring path and tried to make it cheaper with caches, tiles and admission control.
Its own experiments showed that this did not change the cost structure. The frame cache, tiled extraction and scorer micro-optimisations did not give a faster engine.
See the archived [redesign decision](archive/legacy-wp-plan-20261107/redesign-decision-20261103.md).

The rewrite replaces the whole dataflow:

```
raw run (.d / mzML / Thermo .raw)
   │
   ▼
OpenDIACalibrate      identification-free RT, fragment m/z and 1/K0 calibration (no iRT kit, no IDs)
   │  calibration.json, anchors.tsv
   ▼
OpenDIA -calibration  streaming sparse-grid search of the full library inside calibrated windows
   │  candidate peak groups (multi-apex), search features
   ▼
OpenDIACNN crossfit   now (frozen base): trace tensors → CNN logit → stage-2 LDA over CNN + heuristic
                      features, family cross-fitting, pooled paired target-decoy competition → q-values
OpenDIAScoringModel   planned: one voxel model replaces CNN + LDA; same folds and the same paired TDC
```

## Current HeLa result (design run, not independent evidence)

This section records the state of 7–8 November 2026. For the current Astral state, see "Current state" above and [BOARD.md](BOARD.md).

| HeLa 200 ng diaPASEF, paired TDC, experimental q ≤ 0.01 | Targets | SVM-baseline overlap |
|---|---|---|
| OpenSWATH + PyProphet SVM (reference) | 19132 | n/a |
| OpenSWATH + PyProphet XGBoost (reference) | 19725 | n/a |
| Rewrite, no RT prior | 6531 | 5846 |
| Rewrite, CNN-anchor RT window | 8554 | 7502 |
| Rewrite, ID-free calibration filters candidates after extraction | 10679 | 9146 |
| **Rewrite, calibrate then extract (current)** | **11830** | **10158** |
| **Rewrite, teacher-free: panel-40000 calibration anchors + `-ss_iterations 1`** | **13611** | **11533** |
| **Rewrite, teacher-free end to end from one panel-40000 calibration (frozen design)** | **13253** | **11275** |

The current pipeline extracts 198,915 candidates in 40 s and runs CNN crossfit in 1.2 min with a 2.1 GB peak.

**Limits:**
- These are decoy-count estimates on the run used to design the pipeline. They are not a calibrated FDR.
- Every row except the last two uses CNN seeds from the same-run OpenSWATH teacher. The last two use no OpenSWATH input and exceed the teacher-seeded count. Neither has been checked by entrapment.
- On the independent 500 pg run, the CNN alone gave 0 IDs at 1%, paired TDC gave 230, and the calibration failed closed with 0 anchors.
- On the independent ProteoBench Astral entrapment run (raw 01, frozen teacher-free design applied once), the calibration also failed closed with 0 anchors, so no IDs and no entrapment FDP. The HeLa design does not yet transfer to weak or no-1/K0 runs.
- Option B then used raw 01 as a design run, adding a windowed second calibration pass and pooled paired TDC (both opt-in). The frozen result was evaluated once on raw 02 against a hashed gate, and it **failed**.
  - It gave 604 anchors and 360 targets at q ≤ 0.01, with 5 entrapment hits.
  - Combined entrapment FDP estimate 0.028 against a 0.015 limit (95% interval 0.009 to 0.065).
  - On this weak Astral run the experimental q is not conservative enough. See R-WP2 and R-WP4.
- An opt-in second-stage LDA (`OpenDIACNN -stage2`) combines the CNN score with co-elution, library correlation and MS1 evidence.
  - Raw 01 design run: 333 targets at q ≤ 0.01, against 321 for CNN-only.
  - Spent raw 02 re-score (a diagnostic, not a gate): 518 targets, FDP_c 0.012.
  - It has not been evaluated independently. Raw 03 is held back for that.
- The rewrite is about 5,500 to 6,100 precursors short of the reference. Closing that gap is R-WP5.

## Work packages

| WP | Scope | State, 10 November 2026 |
|---|---|---|
| [R-WP1](work-packages/R1-streaming-search.md) | Streaming sparse-grid candidate search | Delivered. Thermo `.raw` streams natively. Opt-in boundary picker. Recall gates open. |
| [R-WP2](work-packages/R2-calibration.md) | ID-free run calibration (RT, m/z, 1/K0, tolerances) | Delivered on HeLa and on Astral raws 01/02 (windowed second pass). Fails on 500 pg. |
| [R-WP3](work-packages/R3-trace-cnn.md) | Raw-trace scorer | CNN + LDA is the frozen base. Work on CNN + LDA is stopped. **Active: single voxel model `OpenDIAScoringModel`.** |
| [R-WP4](work-packages/R4-error-control.md) | Cross-fitting, paired TDC, independent error control | Partial. q ≤ 0.01 holds on raws 01/02. q ≤ 0.05 is anti-conservative. |
| [R-WP5](work-packages/R5-id-gap.md) | Close the gap to OpenSWATH + Percolator | **Active.** Raw 02 at OSW level, raw 01 at 89 %. |
| [R-WP6](work-packages/R6-acceptance.md) | Independent runs, native Thermo, packaging, merge | Open. Raw 03 is held back. |

Dependency order: R-WP5 uses R-WP1 to R-WP4. R-WP6 needs R-WP4 and R-WP5.

**Step board:** [BOARD.md](BOARD.md) lists each step with its state (done, in progress, next, blocked) and the link to its record.

## Rules carried over

- Defaults stay byte-identical. New behaviour is opt-in until it is accepted. Check HeLa default equivalence in a fresh directory before each benchmark run.
- Write the pass rule in a PLAN file (with sha256 and time) before each run or training. Do not change a rule after a result.
- Design on raws 01 and 02 only. Do not tune on an evaluation run. Raw 03 needs approval.
- Entrapment labels never go into training. DIA-NN and OSW are references only. Do not reuse OSW scoring code.
- Never report a calibrated FDR from decoy counts alone. Call q-derived counts decoy-count estimates. Agreement with OpenSWATH is not truth.
- Native Bruker and Thermo input, with no mandatory conversion.
- Inputs and baseline OSW files are read-only.

## Issues

The epic is [#10661](https://github.com/OpenMS/OpenMS/issues/10661). The legacy WP issues #10663 to #10671 and #10675 are closed, each pointing to this plan. #10662 (WP-00) was closed earlier.
See [issue-posting.md](issue-posting.md) for the issue of each R-WP.
