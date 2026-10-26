# WP-05: Exact-query shared extraction and calibration adapter

- **Track:** data handling.
- **Status:** in progress. First read-only extraction increment verified 26 October 2026, full sharing/calibration acceptance still open. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10667](https://github.com/OpenMS/OpenMS/issues/10667).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Share identical extraction queries in the existing engine without changing evidence or waiting for new retrieval science.

## Scope

- Own the adapter to existing mass/RT/IM calibration, with immutable per-pass snapshots and cache invalidation on revision.
- Implement complete query equality, collision-safe sharing and immutable sparse evidence with exact axes and validity masks.
- Retain raw-region access, acquisition geometry and complete picker input/state where existing-engine equivalence requires it.
- Plug extraction opt-in into the existing engine and compare sharing on/off before new candidate/family scoring.
- Handle shared tophat first/last-centroid bugs as a standalone core fix on `develop`, with both IM and non-IM tests, not a compatibility mode.
- Keep that mathematical correction separately reviewed and use the same corrected reference for adapter comparisons.

## First increment: leased read-only extraction

The immediate implementation is **not** a query cache. The full TimSim PQP has 4,000 distinct precursor m/z/RT/IM tuples, so precursor deduplication does not justify that cache on this input. The existing extractor only reads numerical arrays but currently requests a charged, detached mutable copy of every tiled scan.

1. Add a const-array spectrum view to `ISpectrumAccess`, with a default implementation using the existing spectrum method. Existing calibration wrappers therefore keep their transformation and copy behavior.
2. Override that view in private `RawAccess`, aliasing cached arrays through their scan lease. Route `ChromatogramExtractorAlgorithm` through the view. Do not change scoring, extraction mathematics, query identity or scheduling.
3. Extend the existing native raw-access test. A cache that admits a scan but rejects its mutable copy must still extract exactly equal IM/non-IM evidence. Check array identity, last-owner/cancellation lifetime and quadratic-calibration fallback.
4. Run the unchanged full-library tiled 256 MiB Release command once, changing output paths only. Compare the export, all eight feature/score tables and three matrices with the saved normal reference, including finite nonzero IM. Record wall time, RSS, cache counters, private peak and final live bytes.

**Stop rule:** one failing-before/passing-after regression, affected native tests and one successful full workflow. A single timing observation is not a speedup claim. No new dependency, query cache, whole-process memory claim, automatic low-cap scheduling fix, repeated balanced matrix or issue-publication loop in this increment. The remaining WP-05 merge criteria below stay open.

### Observed result, 2026-10-26

Implemented on `feature/opendia-wp05-readonly-extraction`, commit `6423194fc8d8e165003760b13b5f8a737c778762`, stacked on WP-03 `4024bc8325d4fecc35a300891ec7937d10c35624`. The implementation changes six files, adds no dependency or CLI option and leaves mutable scoring access unchanged. The new virtual interface requires a coherent source rebuild, not a binary-compatible drop-in for previously compiled SDK clients. The actual Release tool, libraries and affected tests were rebuilt together.

| Requirement | Concrete check and observed result |
| --- | --- |
| Avoid the unnecessary extraction copy | Added `readOnlyExtraction` to `OpenDIARawAccess_test`. Before the implementation, the actual native extractor fails with `Mutable scan exceeds maximum scan bytes`. Afterwards, the same tight scan/copy cap admits extraction, while a mutable request still fails. Array addresses equal cached column addresses and committed bytes do not increase. |
| Preserve extraction evidence and axes | Tight-budget IM and non-IM extraction matches native reference RT/intensity arrays exactly. The existing real mzML/IM/FAIMS loader comparisons also pass for 1, 16 and 4,096 queries at both cache limits. |
| Keep ownership and cancellation safe | Each of m/z, intensity and mobility independently retains the scan charge after cache cancellation and reader destruction. New admission fails after cancellation. The last array owner releases the charge to zero. |
| Preserve calibrated wrappers | Quadratic and ppm wrapper clones return transformed m/z, not the raw alias. Their native extraction matches direct calibrated readers exactly. Wrappers deliberately retain the existing mutable-copy fallback. |
| Exercise affected native interfaces | `OpenDIARawAccess_test`, `OpenDIARawFailure_test`, `ChromatogramExtractorAlgorithm_test`, `OpenSwathScoring_test` and `SpectrumAccessQuadMZTransforming_test` all pass in the actual Release build. |
| Preserve the full public workflow | One full TimSim tiled256 CLI run finishes with engine/wrapper exit 0. Only output/work paths differ from the previous tiled256 command. Input, library, comparator, command semantics and before/after source/binary/dependency provenance checks pass. |
| Preserve every scientific output | Export: 4,688 rows. FEATURE / MS1 / MS2: 4,700 each. FEATURE_PRECURSOR: 14,219. FEATURE_TRANSITION: 28,200. SCORE_MS2 / PEPTIDE / PROTEIN: 4,688 / 3,287 / 1,871. Precursor / peptide / protein matrices: 985 / 985 / 251. All twelve outputs match the saved normal reference, including IDs, axes, finite nonzero IM, scores and quantities, under the established normalization/tolerance contract. |
| Observe resources without overstating them | Wall 367.38 s, peak RSS 1,039,580 KiB. Private peak 201,325,194 <= 201,326,592 usable bytes, final live 0. Cache hits / misses / evictions: 1,955 / 11,290 / 11,272. This is one observation, not a process-memory bound or balanced comparison. |

The preceding tiled256 observations were 373.04 / 374.40 / 417.98 s and 982,188 / 1,065,304 / 946,476 KiB RSS. The new run proves compatibility and the regression proves copy avoidance, but **does not establish an overall performance or memory improvement**. Its raw misses remain within the preceding range, so repeated decoding is not solved by this change. No further repeats were run.

Local evidence artifacts are retained but not published as a sanitized bundle. `regression-before.log` retains the pre-fix failure, `native-tests-final.log` records the five native checks, and `tiled256/results.json` plus `acceptance.log` record the complete real CLI command, provenance, measurements and each output comparison. `run_readonly.py` reuses the previous strict real-output comparator and provenance helpers rather than a substitute workflow. WP-02 whole-job admission, WP-03 broader acceptance and the remaining WP-05 sharing/calibration criteria remain open.

## Next increment

Follow the [master's ordered next actions](../README.md#current-state-and-next-actions). First attribute the remaining misses to calibration, MS1 extraction, MS2 extraction and scoring in one unchanged full-input run. Reconcile integer event/byte counts with whole-run totals, report overlapping stage times separately, then change only the dominant reread path. Acceptance is one native regression plus one full twelve-output comparison, with no speedup claim from a single observation. Query caching, new retrieval science and full work-package acceptance are not implied by the read-only increment.

## Dependencies

- **Start / local merge:** WP-01 only. No wait on WP-06 for calibration.
- WP-03 tile and WP-04 store adapters integrate as available. Existing inputs support initial extraction work.
- WP-02 hardening is not a blocker. Feed calibrated evidence and raw revisits to WP-06/07.

## Local merge criteria

- [ ] Hand-computed tests cover first/last centroids, empty/singleton scans, duplicate mz and exact/near boundaries for both IM and non-IM tophat extraction.
- [ ] Every query-key field and forced hash collisions are tested, and changed calibration snapshots invalidate affected evidence.
- [ ] Existing-engine sharing on/off runs retain identical IDs and tolerance-equal evidence, scores and quantities against the corrected reference.
- [ ] Repeated consumers share one payload, while lifetime/cancellation tests release buffers only after the last lease.
- [ ] Tile-edge, non-monotone interval, missing-acquisition and FAIMS fixtures match scalar reference extraction without coordinate collapse.
- [ ] The separate core-fix revision and its IM/non-IM test results are recorded independently from the private adapter change.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
