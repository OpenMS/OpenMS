# OpenDIA master plan "Current state", 22 November 2026 (archived 23 November 2026)

Kept unchanged for the record. The current state is in [../../README.md](../../README.md). Relative links are relative to the master plan directory.

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
- **Candidate arm:** `-valley_fraction 0.7 -max_peaks 6` (V7M6). After an output-identical speed-up of the interference resolver (OpenMS `955516d84c`), rule B passed on matched nodes (RESULT_V7M6T2, D r ≈ 7 wall 1.375×). The frozen 5-seed PLAN_V7M6 then **failed on wall only** (RESULT_V7M6_D). D median N @ FDP_c 1 % went from 21,396 to 24,249 (+13.3 %, every seed +11.7 to +16.6 %, 93 % of DIA-NN) with unchanged q bounds, but the median wall was 1.515× against 1.5× on unmatched nodes. The defaults stay; the next step needs a decision.

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
