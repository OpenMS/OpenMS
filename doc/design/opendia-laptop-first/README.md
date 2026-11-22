# OpenDIA rewrite: master plan

**Status, 22 November 2026.** This plan replaced the laptop-first WP-00 to WP-10 plan on 7 November 2026. The old plan is closed and kept unchanged under [archive/legacy-wp-plan-20261107](archive/legacy-wp-plan-20261107/README.md). The status sections up to 20 November are in [archive/status-history-20261122](archive/status-history-20261122/README.md).

The rewrite lives on [`OpenMS/OpenMS:feature/opendia-clean-rewrite`](https://github.com/OpenMS/OpenMS/tree/feature/opendia-clean-rewrite).
All of its code is under `src/topp/OpenDIA/`, apart from the Thermo read-ahead change, which is split out as OpenMS#10925 (approved).
- **Local commits:** the default described below is OpenMS `bcb156e036`, and the dead-code cleanup follows it. Both are local and not yet pushed.
- **Status:** it is an experimental branch, not a pull request. It must not merge into `develop` until the acceptance gates of R-WP6 pass.

The branch's [`src/topp/OpenDIA/README.md`](https://github.com/OpenMS/OpenMS/blob/feature/opendia-clean-rewrite/src/topp/OpenDIA/README.md) describes the tools. The evidence is in [`gates/`](gates/).

## Current state, 22 November 2026

**Benchmark.** HeLa diaPASEF entrapment, r ≈ 7. The library is 1.39 M targets with predicted decoys and homolog-free foreign entrapments. There are two runs:
- **D:** 2737.d, the design run.
- **H:** 22474.d, the held-out run.

The metric is N, precursors at true FDP 1 % (FDP_c), as the median of 5 seeds. DIA-NN 2.2 on the same library gives 26,006 (D) and 41,272 (H).

| Default after | Record | D N @ FDP_c 1 % | H N @ FDP_c 1 % | D / H share of DIA-NN |
|---|---|---:|---:|---|
| Predicted decoys (OpenDIA `-decoys predict`), start of r ≈ 7 | RESULT_ODP_D, RESULT_EXTRA_H (CUR arm) | 8,233 | 20,745 | 32 % / 50 % |
| Evidence block (MS1, mobility, mass error) | RESULT_EXTRA_D/H | 11,332 | 25,643 | 44 % / 62 % |
| Narrow 1/K0 window | RESULT_IMNARROW_D/H | 11,815 | 26,647 | 45 % / 65 % |
| 10 ppm trace tensor | RESULT_TRPPM_D/H | 13,985 | 29,286 | 54 % / 71 % |
| Matched-filter features | RESULT_MT_D/H | 14,895 | 29,873 | 57 % / 72 % |
| Per-run RT refinement | RESULT_RTR_D/H | 15,369 | 30,618 | 59 % / 74 % |
| **Charge-state sibling rescoring** | RESULT_SIB_D/H | 19,704 | 34,871 | 76 % / 84 % |
| **Hard-negative GBM rescoring** | RESULT_GBT_D/H | **21,396** | **37,342** | **82 % / 90 %** |

All records are in `gates/single-binary/pdec/`. Each row passed a frozen 5-seed plan on D. H was run only after the D verdict was committed. Each default flip was byte-identical to the tested arm, and its `-no_*` switch was byte-identical to the previous default. Since the 22 Nov cleanup these stages are fixed: the code has no switch to turn them off.

**What the default does now.**
1. **Calibration:** identification-free calibration.
2. **Search:** a streaming candidate search, 5 peaks per precursor.
3. **Evidence pass:** 6-fragment traces at ±10 ppm, a 12-fragment grid, shared-signal interference, MS1 and mobility evidence, and matched-filter features.
4. **Network scorer:** the `CrossFitScorer` network, with family cross-fitting and self-training.
5. **Sibling rescoring:** four charge-state sibling features, built from out-of-fold first scores, then a rescore. Siblings share the family fold.
6. **GBM rescoring:** three rounds of a cross-fitted histogram GBM, trained on q ≤ 0.01 targets against all windows of the hardest 2 % of decoy precursors.
7. **RT refinement:** a per-run ridge RT correction, then a second full pass.

**Calibration of the tool's own q (pooled at q ≤ 0.01):**

| Run | Lower bound | Combined bound | Pair fairness |
|---|---|---|---|
| D | 1.11 % | 1.27 % | 0.544 |
| H | 1.29 % | 1.47 % | 0.540 |

- **H loosened with the GBM:** it was 1.15 / 1.31 % before (RESULT_GBT_H). The extra IDs are real at true FDP 1 %, but on H the reported q is slightly optimistic. This is the main known cost to watch.
- **Cost:** D takes about 20 min and 35 GB, and H about 30 min and 29 GB, on 8 threads (codon, node-noisy). DIA-NN takes 3-4 min and under 7 GB.

**Where the remaining gap is** (RESULT_GAP2, measured on the sibling default before the GBM):
- **D, scoring:** the DIA-NN-labelled GBM ceiling is now below the tool. 64 % of D's missed IDs at r ≈ 7 have the right peak but score too low.
- **H, candidates:** 38 % of H's missed IDs have no candidate within 6 s.
- **Candidate arm:** `-valley_fraction 0.7 -max_peaks 6` (V7M6) gained D +13.3 % and H +3.6 % at r ≈ 7. It failed rule B only on D r ≈ 7 wall time (1.52× against 1.5×). It needs a user decision: re-time it on matched nodes, or find a cheaper split.

**Removed (22 Nov cleanup).** The code now carries only the best default. Their records stay in `gates/`.

Options that failed their frozen plans or were never adopted:
- `-pg_scores`, `-mobility_apex`, `-faint_features`, `-im_peaks`, `-rt_refine_anchors`, `-hard_decoys` and `-scorer_self_q`;
- the `-trace_tolerance_ppm 0` fixed 0.02 Da window (the value must now be positive);
- the optional LibTorch `OpenDIACNN` path;
- the standalone `OpenDIAInterference` diagnostic tool (the workflow's built-in interference features stay);
- the Python trace-window pilot.

Switches back to an older default, now that each stage has passed:
- `-no_extra_scores`, `-no_matched_filter`, `-no_sibling_features`, `-no_gbm_rescore` and `-gbm_hard_fraction` (fixed at 0.02);
- `OpenDIACalibrate -max_panel_targets` (panel escalation, never used by the workflow);
- code no tool calls: `InterferenceCollector`, `readTraceWindows` and `peakBoundaryName`, found by a `--gc-sections` link of all three tools.

Kept: `-no_rt_refine` (about half the wall time, a real speed choice), the tuning knobs (`-max_peaks`, `-valley_fraction`, `-im_tolerance_scale`, `-trace_tolerance_ppm`), calibration reuse, and the step tools `OpenDIA` and `OpenDIACalibrate`, which the diagnostics use.

The cleanup is OpenMS `319f043f77` (30 files, 4,912 lines removed net) and `9ba060f978` (19 files, 409 lines removed net). 25 of 25 unit tests pass; two CLI tests went with `-max_panel_targets`.

**Astral.** Unchanged since 16 November: OpenDIA reaches 81-84 % of DIA-NN at FDP_c 1 % (archived snapshot). The diaPASEF levers above have not been re-measured on Astral.

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
evidence pass         6-fragment traces (±10 ppm), apex b/y ions, 12-fragment grids, interference graph, MS1 and
                      mobility evidence, matched filter (45 features)
CrossFitScorer        decoy-whitened network, family cross-fitting, self-training
sibling rescoring     4 charge-state sibling features from out-of-fold first scores, network rescored (49 features)
GBM rescoring         3 cross-fitted histogram-GBM rounds against the hardest 2 % of decoy precursors
                      → pooled paired TDC → q-values; then per-run RT refinement and a second full pass
                      (all inside OpenDIAWorkflow)
```

## Work packages

| WP | Scope | State, 22 November 2026 |
|---|---|---|
| [R-WP1](work-packages/R1-streaming-search.md) | Streaming candidate search | Delivered: native-scan XIC peak picking and Thermo `.raw`. **Open: candidate recall on H** (38 % of missed DIA-NN IDs have no candidate; V7M6 awaits a decision). |
| [R-WP2](work-packages/R2-calibration.md) | ID-free run calibration | Delivered, with the windowed pass and per-run RT refinement (default). |
| [R-WP3](work-packages/R3-trace-cnn.md) | Scorer | **C++ `CrossFitScorer` network, then charge-state sibling rescoring, then a hard-negative histogram GBM** (all default). The CNN + LDA path is deleted. |
| [R-WP4](work-packages/R4-error-control.md) | Error control | Predicted decoys. Own q ≤ 0.01 true FDP on r ≈ 7: D 1.11-1.27 %, H 1.29-1.47 % (H loosened with the GBM; watch). |
| [R-WP5](work-packages/R5-id-gap.md) | ID gap | **Active.** HeLa r ≈ 7: **82 % (D) / 90 % (H) of DIA-NN** at FDP_c 1 %. Astral 81-84 % (16 Nov, not re-measured). |
| [R-WP6](work-packages/R6-acceptance.md) | Acceptance, packaging, merge | One-command tool delivered. Dead-code cleanup 22 Nov. Raw 03 held back. PR split open. |

Dependency order: R-WP5 uses R-WP1 to R-WP4. R-WP6 needs R-WP4 and R-WP5.

**Step board:** [BOARD.md](BOARD.md) lists each step with its state (done, in progress, next, blocked) and the link to its record.

## Rules carried over

- Defaults stay byte-identical. New behaviour is opt-in until it is accepted. Check HeLa default equivalence in a fresh directory before each benchmark run.
- Write the pass rule in a PLAN file (with sha256 and time) before each run or training. Do not change a rule after a result.
- Design on raws 01 and 02 only. Do not tune on an evaluation run. Raw 03 needs approval.
- Entrapment labels never go into training. DIA-NN and OSW are references only. Do not reuse OSW scoring code.
- **No OpenSWATH code in OpenDIA.** Use only the OpenMS readers and low-level classes (spectra, chemistry, SQLite).
  - The ctest `OpenDIA_no_openswath` fails on any OpenSWATH reference in the core (`include/`, `source/`).
  - Since 1f29ce2da8 the core is clean: its own PQP reader and library types (`OpenDIALibrary`), its own peak-map
    struct, and the isotope pattern taken directly from `CoarseIsotopePatternGenerator`.
  - Two marked edges remain in tool files. One is `OpenDIAWorkflow -decoys predict` (`OpenSwathLibraryPreparation`),
    which is pending a decision. The other is the `-out_xipm` writer of the `OpenDIA` dev tool, an OpenMS format
    writer that takes OpenSWATH types.
  - OpenDIAWorkflow still links the full libOpenMS (105 MB), so the binary is not yet small.
- Never report a calibrated FDR from decoy counts alone. Call q-derived counts decoy-count estimates. Agreement with OpenSWATH is not truth.
- Native Bruker and Thermo input, with no mandatory conversion.
- Inputs and baseline OSW files are read-only.

## Issues

The epic is [#10661](https://github.com/OpenMS/OpenMS/issues/10661). The legacy WP issues #10663 to #10671 and #10675 are closed, each pointing to this plan. #10662 (WP-00) was closed earlier.
See [issue-posting.md](issue-posting.md) for the issue of each R-WP.
