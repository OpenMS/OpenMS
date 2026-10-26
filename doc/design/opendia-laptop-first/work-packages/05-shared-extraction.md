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

## Native vendor input increment

Commits `8743497d0b` and `3e5852085f` extend the existing adapter without a new dependency, CLI mode or converted-file format. Thermo already had native random access. Bruker now exposes it through the same `IMSExperimentRandomAccess` interface, reusing the existing decoders and default streaming calibration/processing. The tiled CLI no longer rejects vendor input before dispatch. Both readers are enabled in the actual Release build and remain required for validation/distribution. No internal indexed-mzML conversion is required or performed.

| Requirement | Actual check and result |
| --- | --- |
| Build and public native interfaces | Coherent Release build with `WITH_OPENTIMS=ON` and `WITH_THERMO_RAW=ON`. All six focused CTests pass: raw access/failure, Bruker, Thermo random access, SWATH loading and on-disc access. Optional downloaded HeLa fixture sections are not claimed as executed. |
| Bruker native semantics | Real frozen TimSim `.d.zip`, 25 maps and 75 sampled first/middle/last scans. Exact isolation/IM geometry, map/scan membership, RT/MS level and all numeric arrays match the independent current streaming reader. Cached views and sampled IM/non-IM extraction match native access; final private live/reserved bytes are zero. The additional all-spectrum streaming parity regression is retained behind the existing optional DIA fixture guard. |
| Thermo native semantics | Real 2.9 GB Astral `.raw`, 301 maps and 903 sampled scans. Native metadata, arrays, cached views, RT lookup and sampled non-IM extraction match; final private live/reserved bytes zero. This is not a complete Thermo identification/scoring acceptance run. |
| Complete native Bruker workflow | One complete tiled256 run and one same-input native normal reference. Same binary, full 4,000-precursor / 24,000-transition library, scientific options and threads. All twelve outputs match under the existing ID/tolerance rules, with finite nonzero IM. Export 4,699. FEATURE / MS1 / MS2 4,716 each, PRECURSOR 14,409, TRANSITION 28,296. SCORE_MS2 / PEPTIDE / PROTEIN 4,699 / 3,287 / 1,871. Precursor / peptide / protein matrices 942 / 942 / 251. |
| Private accounting versus process resources | Tiled: 327.32 s, 966,944 KiB peak RSS, private peak 201,326,507 <= 201,326,592 usable bytes, final live zero. Hits / misses / evictions 1,844 / 11,433 / 11,415. Native normal: 186.31 s, 1,522,012 KiB RSS. These are single observations, not a balanced speedup or whole-process bound. Thermo's reader check alone reached 4,352,800 KiB RSS in 25.06 s, so native memory acceptance remains open. |
| No hidden conversion or discarded failure | Both real vendor checks open the original vendor path. The full native workflow calls no converter. The older converted-mzML comparison failed, 4,688 versus 4,699 exported rows and 985 versus 942 precursor/peptide matrix rows. It remains a failed cross-representation comparison. It is not explained away by tolerance or claimed to match. |

Bruker startup decodes frames once to preserve the streaming reader's exclusion of empty MS2 windows and its MS1 metadata. It retains metadata/index membership, not a converted file. Archive ownership survives clones. A shared native handle serializes clone decoding, an explicit throughput ceiling rather than a parallel-reader claim. Native decoding/index/workspace memory is outside the private numeric cap.

The first clean-service Thermo check failed because the pre-existing nonstandard .NET installation was not discoverable. Supplying its `DOTNET_ROOT` made the real native check pass. A distribution must install/discover its supported runtime without asking users to configure this or convert their data. Clean-machine runtime packaging, full Thermo science, the precise cross-representation difference and broader native memory/acquisition gates remain open.

Local evidence is in `wp05-native-20261026`: `native-tests.log`, `streaming-check.log`, `native-thermo.log`, and the two `native-bruker-*/results.json` records with command/source/binary/input/dependency provenance. The original tiled record deliberately retains its failed converted-reference status. `native-reference-acceptance.log` records the successful twelve-output native normal/tiled comparison. No tiled repeat, benchmark matrix or package closure followed.

## Next increment

Follow the [master's ordered next actions](../README.md#current-state-and-next-actions). Keep both native readers enabled and use the native-input reference. First attribute the remaining misses to calibration, MS1 extraction, MS2 extraction and scoring in one unchanged full-input native Bruker run. Separately, one bounded Thermo reader check should record live native handle counts and RSS after indexing, map construction, raw-adapter construction and release before changing handle ownership. Reconcile integer event/byte counts with whole-run totals, report overlapping stage times separately, then change only the dominant reread path. Acceptance is one native regression plus one full twelve-output comparison, with no speedup claim from a single observation. Query caching, new retrieval science and full work-package acceptance are not implied by the read-only increment.

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
