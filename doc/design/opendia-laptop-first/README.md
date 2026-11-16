# OpenDIA rewrite: master plan

**Status, 16 November 2026.** This plan replaced the laptop-first WP-00 to WP-10 plan on 7 November 2026. The old plan is closed and kept unchanged under [archive/legacy-wp-plan-20261107](archive/legacy-wp-plan-20261107/README.md).

The rewrite lives on [`OpenMS/OpenMS:feature/opendia-clean-rewrite`](https://github.com/OpenMS/OpenMS/tree/feature/opendia-clean-rewrite).
All of its code is under `src/topp/OpenDIA/`, and it makes no changes to the shared OpenMS API. The branch tip is `bc30ae9394`.
It is an experimental branch, not a pull request. It must not merge into `develop` until the acceptance gates of R-WP6 pass.
The branch's [`src/topp/OpenDIA/README.md`](https://github.com/OpenMS/OpenMS/blob/feature/opendia-clean-rewrite/src/topp/OpenDIA/README.md) describes the tools. The evidence is in [`gates/`](gates/).

## Current state, 16 November 2026

### Snapshot against DIA-NN (16 Nov, `gates/single-binary/STATUS_vs_DIANN.txt`)
- **Setup:**
  - Current `OpenDIAWorkflow`: self-calibrating, arm-E search library, `library12_full.pqp` scoring library, seeds 0-4,
    codon on 8 threads. Code `ab46f97ad4`.
  - DIA-NN 2.2.0 on the same 302,480 targets with its own decoys (`--gen-decoys`).
  - Both cut at their own run q ≤ 0.01. Entrapment FDP_c is computed only for evaluation.

| | raw 01 | raw 02 | raw 03 (held out) |
|---|---|---|---|
| OpenDIA N at FDP_c 1 % (median, range) | 1804 (1728-1827) | 1909 (1867-1935) | 1907 (1831-1957) |
| DIA-NN N at FDP_c 1 % | 2067 | 2283 | not run |
| OpenDIA / DIA-NN | 87 % | 84 % | |
| OpenDIA IDs at q ≤ 0.01 (median), FDP_c | 1709, 0.0047-0.0095 | 1811, 0.0012-0.0088 | 1849, 0.0058-0.0096 |
| DIA-NN IDs at q ≤ 0.01, FDP_c | 2361, **0.0144** | 2322, **0.0112** | |
| DIA-NN's q ≤ 0.01 IDs that OpenDIA also reports | 67 % (1588) | 72 % (1676) | |
| Same, human (non-entrapment) IDs only | 68 % | 73 % | |
| OpenDIA's q ≤ 0.01 IDs that DIA-NN also reports | 93 % | 93 % | |
| OpenDIA only / DIA-NN only (q ≤ 0.01) | 121 / 773 | 128 / 646 | |
| Overlap at each tool's own FDP_c 1 % cut | 75 % of DIA-NN's | 76 % of DIA-NN's | |

- **At equal true error (FDP_c 1 %)**, OpenDIA reaches 84-87 % of DIA-NN.
- **At each tool's own q ≤ 0.01**, the gap looks larger: 72-78 %. One reason is that DIA-NN's q ≤ 0.01 is
  anti-conservative on this benchmark (FDP_c 1.1-1.4 %), while OpenDIA's is conservative (≤ 1 %).
- **93 % of what OpenDIA reports is also in DIA-NN's list.** The disagreement is mostly DIA-NN IDs that OpenDIA does
  not reach.
- **Where the missed DIA-NN IDs are lost** (current funnel, `gates/extraction-parity/gap/FUNNEL_workflow_library12full.txt`,
  medians, raw 01 / 02):

  | stage | raw 01 | raw 02 |
  |---|---|---|
  | not scored | 8 | 10 |
  | wrong peak | 116 | 106 |
  | right peak, below the cut | 570 | 463 |

  So about 80 % is ranking.

**Run time and resources, one Astral 15-min raw (3.1 GB .raw).** OpenDIA medians come from the 15 codon runs (one per
seed, 3 raws).

| OpenDIA stage | seconds |
|---|---|
| calibration (2 passes) | 128-153 |
| search | 15-19 |
| one evidence pass | 148-240 |
| features | 7-9 |
| scoring | 17-25 |
| **total wall** | **6:08-7:35** |

| | OpenDIA | DIA-NN 2.2 (local docker, same library) |
|---|---|---|
| Wall time | 6:08-7:35 (8 threads) | 13:27 raw 01 (14 threads); 8:48 raw 02 (9 threads, run beside another job) |
| CPU time | 586-714 s | not recorded |
| Peak memory | **21.9-22.7 GB** | not recorded (the docker run hid it) |
| Input | Thermo .raw read directly | mzML (FileConverter, about 5.8 GB per raw), conversion time not included |

- OpenDIA is about 1.3-2× faster than DIA-NN here.
- OpenDIA's memory rose from 9.6 GB (with `library12_pop`) to about 22 GB with the full 12-fragment scoring library.
  That is an open efficiency item.

**Goal.** OpenDIA must identify as many precursors as DIA-NN on the same library, with honest error control and no teacher. The earlier goal was the OpenSWATH + Percolator level, and it is passed.

**Design data.** ProteoBench Astral entrapment raws 01 and 02 are the design runs. Raw 03 is the held-out run: it is evaluated only after a plan's verdict on raws 01/02 is committed (approved by the user on 16 November 2026).

**Metric.** N is the number of target precursors at an entrapment FDP_c of 1 %, with FDP_c = N_e (1 + 1/r) / N. Every count below is a median over five scorer seeds. The seed range is 50-200 IDs.

**One binary: `OpenDIAWorkflow`** (OpenMS `794ee1565e`, README `e612e34d9c`). It runs on Thermo .raw directly, and there is no Python in the analysis. The steps are:
- calibration;
- candidate search;
- one evidence pass, with 24 features (search, boundary, b/y ions, 12-fragment evidence, fragment deficit and shared-signal interference);
- a cross-fitted C++ scorer (a small network with a decoy-whitened readout and self-training);
- pooled paired TDC q-values per run, and a global q-value over runs.

Cost: about 6-7.5 min and 22 GB per Astral raw on 8 threads with the full 12-fragment scoring library (9.6 GB with `library12_pop`). Output is byte-identical across repeats and thread counts.

| N at entrapment FDP_c 1 % (median of 5 seeds) | raw 01 | raw 02 | FDP_c at q ≤ 0.01 |
|---|---|---|---|
| OpenDIAWorkflow, reused calibration | 1819 | 1870 | 0.0024-0.0082 |
| OpenDIAWorkflow, self-calibrating, `library12_pop` scoring library | 1718 | 1941 | raw 01 0.0048-0.0104, raw 02 0.0023-0.0045 |
| OpenDIAWorkflow, self-calibrating, full 12-fragment scoring library `library12_full.pqp` (current setup; held-out raw 03: 1907) | 1804 | 1909 | 0.0012-0.0096 on all 15 runs (raws 01-03) |
| Python chain it replaces (D11FDI, C++ interference features) | 1848 | 1925 | ≤ 0.0071 |
| DIA-NN 2.2, same targets, its own decoys (reference) | 2067 | 2301 | |
| Former native base (CNN + LDA) | 1029 | 1298 | |

- **Decoys.** Since PLAN_DECOYGEN, the decoys come from OpenSwathDecoyGenerator `-predict_decoys`, which gives each decoy its own predicted spectrum and RT. This removed the entrapment excess that failed every earlier arm.
- **Interference features** (PLAN_GAP arm I). These gave +118 to +185 on raw 01 and +199 to +351 on raw 02 at 1 % over seeds. About a third of the raw 02 gain comes from entrapment pairs that only exist in this benchmark.
- **Self-calibrated RT window.**
  - On raw 01 it is wide (151 s).
  - A refined spread fit gets the spread right, but over five seeds it scores +5.8 % on raw 01 and -3.1 % on raw 02, so it fails (PLAN_CALREFINE2).
  - Two label-free width rules are undefined or point the wrong way.
  - A cross-fitted per-run width chosen on IDs (PLAN_RTWIN2) scores −2.2 % on raw 01 and −2.3 % on held-out raw 03, so it fails. That post-hoc filter was removed.
  - Recalibrating from the run's own IDs, with a window that follows the local RT error (PLAN_RTPROFILE), recovers the gradient ends: +23 to +24 IDs early, +7 to +17 late. But the narrower middle loses (raw 02 −52), the strict-FDR ranking is less stable, and it fails: raw 02 −1.05 %, raw 01 seed 4 FDP_c 0.0104, held-out raw 03 −0.6 %. Removed (`gates/single-binary/rtwin/`).
  - With the full 12-fragment scoring library, the wide window keeps FDP_c ≤ 0.0096 on all 15 runs, so the width is no longer an FDR problem.
- **Gap to DIA-NN.** It is 13-16 % at FDP_c 1 %, almost entirely in score formation on the right peak (`gates/extraction-parity/gap/`).
  - A joint (presence-weighted) interference fit has no room. Block I's explaining-away comes from claimants that are themselves absent (PRECHECK_JOINT).
  - A labelled ceiling on today's 24 inputs recovers only about +100 per raw at today's decoy rate.
  - Scorer-capacity attempts on those inputs all failed:
    - signal-level interactions and self-trained trees, offline;
    - a 3-seed ensemble (PLAN_ENSEMBLE): +0.2 % / +0.6 %, raw 03 +0.3 %. It does cut the raw 01 seed range from 99 to 17.
  - Closing the gap needs new evidence for the weak right-peak IDs (next).

**Direction.** The single voxel model was replaced in practice by the heuristic-feature network (D11 line, R-WP3). That network is now in C++ as `OpenDIA::CrossFitScorer`. OpenDIACNN and LibTorch are no longer on the default path.

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
evidence pass         6-fragment traces, apex b/y ions, 12-fragment grids, shared-signal interference graph (24 features)
CrossFitScorer        decoy-whitened network, family cross-fitting, self-training, pooled paired TDC → q-values
                      (all inside OpenDIAWorkflow; OpenDIACNN is the earlier, optional CNN + LDA scorer)
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

| WP | Scope | State, 16 November 2026 |
|---|---|---|
| [R-WP1](work-packages/R1-streaming-search.md) | Streaming candidate search | Delivered: native-scan XIC peak picking and Thermo `.raw`. Recall at calibrated windows is open. |
| [R-WP2](work-packages/R2-calibration.md) | ID-free run calibration | Delivered, with the windowed pass. **Open: the self-calibrated RT window width** (raw 01 too wide). |
| [R-WP3](work-packages/R3-trace-cnn.md) | Scorer | **C++ `CrossFitScorer` in `OpenDIAWorkflow`**: 98.6 / 99.0 % of the Python D11FDC median. CNN + LDA stopped. |
| [R-WP4](work-packages/R4-error-control.md) | Error control | q ≤ 0.01 holds with reused calibration (FDP_c ≤ 0.0082). Self-calibrated raw 01 reaches 0.0104 on 2 of 5 seeds. |
| [R-WP5](work-packages/R5-id-gap.md) | ID gap | **Active.** Past OSW + Percolator. 88 / 81 % of DIA-NN (raw 01 / 02, reused calibration). |
| [R-WP6](work-packages/R6-acceptance.md) | Acceptance, packaging, merge | One-command tool delivered. Raw 03 held back. PR split open. |

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
