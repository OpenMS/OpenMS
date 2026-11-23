# OpenDIA rewrite: master plan

**Status, 23 November 2026.** This plan replaced the laptop-first WP-00 to WP-10 plan on 7 November 2026. The old plan is closed and kept unchanged under [archive/legacy-wp-plan-20261107](archive/legacy-wp-plan-20261107/README.md). Earlier status sections are in [archive/status-history-20261122](archive/status-history-20261122/README.md) (to 20 November) and [README_state_20261122.md](archive/status-history-20261122/README_state_20261122.md) (22 November).

The rewrite lives on [`OpenMS/OpenMS:feature/opendia-clean-rewrite`](https://github.com/OpenMS/OpenMS/tree/feature/opendia-clean-rewrite).
All of its code is under `src/topp/OpenDIA/`, apart from the Thermo read-ahead change, which is split out as OpenMS#10925 (approved).
- **Local commits:** the default described below is OpenMS `bcb156e036`. The cleanup (`319f043f77`, `9ba060f978`, `35db7eaee1`), the faster interference resolver (`955516d84c`) and the rename to `OpenDIA` (`46243afae4`) follow it, all output-identical. The GBM final-round fix (`88e5da3023`) follows too; its HeLa byte-check is pending. None of them is pushed yet.
- **Status:** it is an experimental branch, not a pull request. It must not merge into `develop` until the acceptance gates of R-WP6 pass.

The branch's [`src/topp/OpenDIA/README.md`](https://github.com/OpenMS/OpenMS/blob/feature/opendia-clean-rewrite/src/topp/OpenDIA/README.md) describes the tool. The evidence is in [`gates/`](gates/).

## Current state, 23 November 2026

**Benchmark.** HeLa diaPASEF entrapment, r ≈ 7. The library is 1.39 M targets with predicted decoys and homolog-free foreign entrapments. There are two runs:
- **D:** 2737.d, the design run.
- **H:** 22474.d, the held-out run.

The metric is N, precursors at true FDP 1 % (FDP_c), as the median of 5 seeds. DIA-NN 2.2 on the same library gives 26,006 (D) and 41,272 (H).

| Default after | Record | D N @ FDP_c 1 % | H N @ FDP_c 1 % | D / H share of DIA-NN |
|---|---|---:|---:|---|
| Predicted decoys (`-decoys predict`), start of r ≈ 7 | RESULT_ODP_D, RESULT_EXTRA_H (CUR arm) | 8,233 | 20,745 | 32 % / 50 % |
| Evidence block (MS1, mobility, mass error) | RESULT_EXTRA_D/H | 11,332 | 25,643 | 44 % / 62 % |
| Narrow 1/K0 window | RESULT_IMNARROW_D/H | 11,815 | 26,647 | 45 % / 65 % |
| 10 ppm trace tensor | RESULT_TRPPM_D/H | 13,985 | 29,286 | 54 % / 71 % |
| Matched-filter features | RESULT_MT_D/H | 14,895 | 29,873 | 57 % / 72 % |
| Per-run RT refinement | RESULT_RTR_D/H | 15,369 | 30,618 | 59 % / 74 % |
| **Charge-state sibling rescoring** | RESULT_SIB_D/H | 19,704 | 34,871 | 76 % / 84 % |
| **Hard-negative GBM rescoring** | RESULT_GBT_D/H | **21,396** | **37,342** | **82 % / 90 %** |

All records are in `gates/single-binary/pdec/`. Each row passed a frozen 5-seed plan on D, and H was run only after the D verdict was committed. Each default flip was byte-identical to the tested arm. The code now carries only this default, so the stages cannot be switched off.

**What the default does.** `OpenDIA` is the only tool. Per run:
1. **Calibration:** identification-free calibration of RT, fragment m/z and 1/K0, with a windowed second pass if needed.
2. **Search:** a streaming candidate search, 5 peaks per precursor.
3. **Evidence pass:** 6-fragment traces at ±10 ppm, a 12-fragment grid, shared-signal interference, MS1 and mobility evidence, and matched-filter features.
4. **Network scorer:** the `CrossFitScorer` network, with family cross-fitting and self-training.
5. **Sibling rescoring:** four charge-state sibling features, built from out-of-fold first scores, then a rescore.
6. **GBM rescoring:** three rounds of a cross-fitted histogram GBM, trained on q ≤ 0.01 targets against all windows of the hardest 2 % of decoy precursors.
7. **RT refinement:** a per-run ridge RT correction, then a second full pass.

**Calibration of the tool's own q (pooled at q ≤ 0.01):**

| Run | Lower bound | Combined bound | Pair fairness |
|---|---|---|---|
| D | 1.11 % | 1.27 % | 0.544 |
| H | 1.29 % | 1.47 % | 0.540 |

- **H loosened with the GBM:** it was 1.15 / 1.31 % before (RESULT_GBT_H). The extra IDs are real at true FDP 1 %, but on H the reported q is slightly optimistic. This is the main known cost to watch.
- **Cost:** D takes about 19 min and 34 GB, and H about 30 min and 29 GB, on 8 threads (codon). DIA-NN takes 3-4 min and under 7 GB. This is now the largest practical gap.

**Where the remaining gap is** (RESULT_GAP2, measured before the GBM):
- **D, scoring:** 64 % of D's missed IDs have the right peak but score too low. A DIA-NN-labelled GBM ceiling is below the tool.
- **H, candidates:** 38 % of H's missed IDs have no candidate within 6 s.
- **Candidate arm V7M6** (`-valley_fraction 0.7 -max_peaks 6`, opt-in):
  - The faster interference resolver (`955516d84c`, output-identical; 1,562 → 878 s on the D r ≈ 7 first pass) brought it under the wall bar on matched nodes (RESULT_V7M6T2, 1.375×).
  - The frozen 5-seed PLAN_V7M6 then **failed on wall only** (RESULT_V7M6_D). D median N went from 21,396 to 24,249 (+13.3 %, every seed +11.7 to +16.6 %, 93 % of DIA-NN), with the same q bounds and fairness 0.540. But the median wall was 1.515× against 1.5×, on unmatched nodes.
  - The defaults stay. The next step needs a user decision: keep it opt-in, speed it up and freeze a new plan, or change the timing rule for future plans.

**Code cleanup (22-23 Nov).** The code carries only the best default; the records of everything removed stay in `gates/`.
1. **Failed or unadopted options** (`319f043f77`, 4,912 lines removed net):
   - `-pg_scores`, `-mobility_apex`, `-faint_features`, `-im_peaks`, `-rt_refine_anchors`, `-hard_decoys`, `-scorer_self_q`;
   - `-trace_tolerance_ppm 0`, the LibTorch `OpenDIACNN` path, the `OpenDIAInterference` tool and the Python pilot.
2. **Switches back to older defaults** (`9ba060f978`, 409 lines):
   - `-no_extra_scores`, `-no_matched_filter`, `-no_sibling_features`, `-no_gbm_rescore`, `-gbm_hard_fraction`;
   - `OpenDIACalibrate -max_panel_targets`, plus linker-dead code.
3. **One tool** (`35db7eaee1`, 2,202 lines removed net):
   - The old step tools `OpenDIA` and `OpenDIACalibrate` are gone.
   - Also removed: code only the step tools used (`-out_xipm` peak-map export, the search peak and feature tables, the trace dump, the calibration diagnostics and anchors TSV), the workflow's `-calibration`/`-anchors` reuse, and the old `TOPP_OpenDIA_*` tests of the pre-rewrite front-end (they no longer ran).
   - Earlier diagnostics (PRECHECK_NOCAND, PRECHECK_IDGAP, ...) that used the step tools stay reproducible from their recorded containers.
4. **Rename** (`46243afae4`): the one tool `OpenDIAWorkflow` is now `OpenDIA` (directory `src/topp/OpenDIA/`, binary `OpenDIA`). `-scoring_lib`, removed in step 3, is restored, because the Astral benchmark scores with `library12_full.pqp`. Byte-identical on D r ≈ 1.56 (3/3 files, codon `rename-20261123`).

Kept: `-no_rt_refine` (about half the wall time), the tuning knobs (`-max_peaks`, `-valley_fraction`, `-im_tolerance_scale`, `-trace_tolerance_ppm`, `-trace_tensor_gb`, `-rt_refine_min_ids`), `-scoring_lib` and `-decoys predict`. 14 of 14 OpenDIA tests pass with no build warnings. The one-tool container (`opendia-pd10-onebin.sif`, sha256 `ecf60560…`) is byte-identical to the reference outputs on all 9 files of D r ≈ 1.56, D r ≈ 7 V7M6 and H r ≈ 7 (codon `onebin-20261123`).

**Astral, re-measured 23 Nov on the laptop profile (16 threads, 32 GB; [RESULT_ASTRAL_LAPTOP](gates/single-binary/astral-laptop/RESULT_ASTRAL_LAPTOP.md), interim).**
- **FAIL, robustness:** 5 of 15 runs crashed in the GBM rescore (`a fold has zero decoy spread`, after a degenerate round 3). The fix (`88e5da3023`) falls back to the round 2 scores. The rerun and the HeLa byte-check are pending.
- **The 10 runs that finished:** median N @ FDP_c 1 % is 2041 / 2128 / 2252 (raws 01 / 02 / 03), against 1804 / 1909 / 1907 on 16 Nov. That is 93 / 90 / 100 % of DIA-NN (.raw), up from 81-84 %. But these are survivors only.
- **New concern:** own q ≤ 0.01 FDP_c is above 1 % on 6 of 10 runs (up to 1.77 %). On 16 Nov it was never above 1 %. There are also two low seeds (raw02_s4 1301).
- **Cost:** 6-8 min median wall and about 17 GB peak RSS per raw.

## Why the direction changed

The legacy plan kept OpenSWATH's chromatogram and feature-scoring path and tried to make it cheaper with caches, tiles and admission control.
Its own experiments showed that this did not change the cost structure. The frame cache, tiled extraction and scorer micro-optimisations did not give a faster engine.
See the archived [redesign decision](archive/legacy-wp-plan-20261107/redesign-decision-20261103.md).

The rewrite replaces the whole dataflow:

```
raw run (.d / mzML / Thermo .raw)
   │
   ▼
calibration           identification-free RT, fragment m/z and 1/K0 calibration (no iRT kit, no IDs)
   │  calibration.json
   ▼
search                streaming sparse-grid search of the full library inside calibrated windows
   │  candidate peak groups (multi-apex), search features
   ▼
evidence pass         6-fragment traces (±10 ppm), apex b/y ions, 12-fragment grids, interference graph, MS1 and
                      mobility evidence, matched filter (45 features)
CrossFitScorer        decoy-whitened network, family cross-fitting, self-training
sibling rescoring     4 charge-state sibling features from out-of-fold first scores, network rescored (49 features)
GBM rescoring         3 cross-fitted histogram-GBM rounds against the hardest 2 % of decoy precursors
                      → pooled paired TDC → q-values; then per-run RT refinement and a second full pass
                      (all steps inside the one tool, OpenDIA)
```

## Work packages

| WP | Scope | State, 23 November 2026 |
|---|---|---|
| [R-WP1](work-packages/R1-streaming-search.md) | Streaming candidate search | Delivered: native-scan XIC peak picking and Thermo `.raw`. **Open: candidate recall on H** (38 % of missed DIA-NN IDs have no candidate; V7M6 failed PLAN_V7M6 on wall only and awaits a decision). |
| [R-WP2](work-packages/R2-calibration.md) | ID-free run calibration | Delivered, with the windowed pass and per-run RT refinement (default). |
| [R-WP3](work-packages/R3-trace-cnn.md) | Scorer | **C++ `CrossFitScorer` network, then charge-state sibling rescoring, then a hard-negative histogram GBM** (all default). The CNN + LDA path is deleted. |
| [R-WP4](work-packages/R4-error-control.md) | Error control | Predicted decoys. Own q ≤ 0.01 true FDP on r ≈ 7: D 1.11-1.27 %, H 1.29-1.47 % (H loosened with the GBM; watch). |
| [R-WP5](work-packages/R5-id-gap.md) | ID gap | **Active.** HeLa r ≈ 7: **82 % (D) / 90 % (H) of DIA-NN** at FDP_c 1 %. Astral (23 Nov, interim, 10 of 15 runs): 90-100 %, but own q is loose and 5 runs crashed (fix pending). |
| [R-WP6](work-packages/R6-acceptance.md) | Acceptance, packaging, merge | One tool, `OpenDIA` (renamed 23 Nov); cleanup 22-23 Nov (step tools removed). Raw 03 held back. PR split open. Speed and memory against DIA-NN open. |

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
  - One marked edge remains in the tool file: `OpenDIA -decoys predict` (`OpenSwathLibraryPreparation`),
    which is pending a decision.
  - OpenDIA still links the full libOpenMS (105 MB), so the binary is not yet small.
- Never report a calibrated FDR from decoy counts alone. Call q-derived counts decoy-count estimates. Agreement with OpenSWATH is not truth.
- Native Bruker and Thermo input, with no mandatory conversion.
- Inputs and baseline OSW files are read-only.

## Issues

The epic is [#10661](https://github.com/OpenMS/OpenMS/issues/10661). The legacy WP issues #10663 to #10671 and #10675 are closed, each pointing to this plan. #10662 (WP-00) was closed earlier.
See [issue-posting.md](issue-posting.md) for the issue of each R-WP.
