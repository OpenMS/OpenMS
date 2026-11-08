# OpenDIA rewrite: master plan

**Status, 7 November 2026.** This replaces the laptop-first WP-00 to WP-10 plan, which is now closed and kept unchanged under [archive/legacy-wp-plan-20261107](archive/legacy-wp-plan-20261107/README.md).

The rewrite lives on [`OpenMS/OpenMS:feature/opendia-clean-rewrite`](https://github.com/OpenMS/OpenMS/tree/feature/opendia-clean-rewrite).
That is 30 commits on `develop`, all under `src/topp/OpenDIA/`, with no shared OpenMS API changes.
It is an experimental branch, not a pull request, and it must not merge into `develop` until the acceptance gates of R-WP6 pass.
Detailed evidence lives in the branch's [`src/topp/OpenDIA/README_TRACE_MODEL.md`](https://github.com/OpenMS/OpenMS/blob/feature/opendia-clean-rewrite/src/topp/OpenDIA/README_TRACE_MODEL.md).

## Why the direction changed

The legacy plan kept OpenSWATH's chromatogram and feature-scoring path and tried to make it cheaper with caches, tiles and admission control.
Its own experiments showed that this did not change the cost structure: the frame cache, tiled extraction and scorer micro-optimisations did not give a faster engine.
See the archived [redesign decision](archive/legacy-wp-plan-20261107/redesign-decision-20261103.md).

The rewrite replaces the whole dataflow:

```
raw run (.d / mzML)
   │
   ▼
OpenDIACalibrate      identification-free RT, fragment m/z and 1/K0 calibration (no iRT kit, no IDs)
   │  calibration.json
   ▼
OpenDIA -calibration  streaming sparse-grid search of the full library inside calibrated windows
   │  candidate peak groups (multi-apex)
   ▼
OpenDIACNN crossfit   raw XIC/XIM trace tensors → CNN, family cross-fitting,
                      paired target-decoy competition → precursor q-values
```

## Current HeLa result (design run, not independent evidence)

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

| WP | Scope | State |
|---|---|---|
| [R-WP1](work-packages/R1-streaming-search.md) | Streaming sparse-grid candidate search | Delivered, open for recall gates |
| [R-WP2](work-packages/R2-calibration.md) | ID-free run calibration (RT, m/z, 1/K0, tolerances) | Delivered on HeLa; fails on 500 pg |
| [R-WP3](work-packages/R3-trace-cnn.md) | Raw-trace CNN scorer and native training | Delivered |
| [R-WP4](work-packages/R4-error-control.md) | Cross-fitting, paired TDC, independent error control | Partial: same-run teacher, no entrapment |
| [R-WP5](work-packages/R5-id-gap.md) | Close the gap to OpenSWATH + PyProphet | **Active** |
| [R-WP6](work-packages/R6-acceptance.md) | Independent runs, native Thermo, packaging, merge | Open |

Dependency order: R-WP5 uses R-WP1 to R-WP4. R-WP6 needs R-WP4 and R-WP5.

## Rules carried over

- Defaults stay byte-identical; new behaviour is opt-in until accepted.
- Fix a gate rule before running. Design on HeLa, then apply once and unchanged to an independent run, with no tuning on evaluation runs.
- Never report a calibrated FDR from same-run teacher-seeded estimates. Agreement with OpenSWATH is not truth.
- Native Bruker and Thermo input, with no mandatory conversion (carried over from the legacy plan).
- Inputs and baseline OSW files are read-only.

## Issues

The epic is [#10661](https://github.com/OpenMS/OpenMS/issues/10661). The legacy WP issues #10663 to #10671 and #10675 are closed, each pointing to this plan; #10662 (WP-00) was already closed.
See [issue-posting.md](issue-posting.md) for the new mapping.
