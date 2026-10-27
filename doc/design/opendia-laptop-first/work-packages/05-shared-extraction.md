# WP-05: Exact-query shared extraction and calibration adapter

- **Track:** data handling.
- **Status:** in progress. First read-only extraction increment verified 26 October 2026, full sharing/calibration acceptance still open. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10667](https://github.com/OpenMS/OpenMS/issues/10667).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Provide efficient calibrated evidence access and extraction, including an OpenDIA-owned replacement dataflow where warranted. Existing-engine query sharing is a compatibility experiment, not the final architecture. Follow the [native-region/fused execution experiment](../performance-strategy.md#immediate-experiment-native-regions-through-real-scoring).

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

## Phase-attribution increment

Commit `357128a7b2` adds fixed-size counters to the existing cache and scoped read purposes at actual extraction-worker and scoring boundaries. CLI calibration joins override the purpose during calibration. No new dependency, profiling framework, CLI option or scientific/scheduling change was added. Both native readers remain enabled.

| Requirement | Concrete check and observed result |
| --- | --- |
| Correct phase ownership and accounting | The existing native raw-access test now checks nested/exception/thread isolation, extraction through the quadratic-calibration fallback, MS1/MS2 separation, calibration override, numeric-byte definitions and pressure-eviction ownership. Independent totals and every phase field reconcile. All five focused Release CTests pass. |
| Actual unchanged public workflow | One full native Bruker tiled256 command, 4,000 precursors / 24,000 transitions, eight threads, unchanged input/library/science/cap. Engine and wrapper exit 0, unchanged before/after provenance. Export, eight feature/score tables and three matrices match the saved native normal reference under the existing ID/tolerance contract, including finite nonzero IM. Row counts match the native evidence above. |
| Complete attribution, including failure detection | All ten emitted count/byte/time fields sum exactly to independent totals. Hits/misses/evictions also match the original cache line. All four named phases are exercised, unclassified fields are zero, and a deliberately corrupted total is rejected by the acceptance check. |

| Phase | Decodes | Decoded numeric bytes | Mutable copies | Decoder time, seconds |
| --- | ---: | ---: | ---: | ---: |
| Calibration | 2,291 | 5,499,196,440 | 1,163 | 23.066 |
| MS1 extraction | 418 | 3,523,962,408 | 744 | 8.471 |
| MS2 extraction | 703 | 1,134,799,416 | 703 | 7.442 |
| Scoring | 7,978 | 57,461,919,792 | 9,432 | 178.332 |
| Total | 11,390 | 67,619,878,056 | 12,042 | 217.312 |

Scoring owns **70.04% of decodes, 84.98% of decoded numeric bytes and 82.06% of decoder time**. This measures ownership, not a safe duplicate key or an MS1/MS2 split within scoring. Numeric bytes count materialized double payloads, not compressed/native frame/index/workspace bytes. Per-call read/decode/copy durations overlap across threads; decode time is inside read time. Total read time 340.896 s is not exclusive stage wall time and must not be added to decoder time.

Observed wall 291.96 s, peak RSS 1,033,472 KiB, private peak 201,326,215 <= 201,326,592 usable bytes, final live zero. Hits / misses / evictions: 1,887 / 11,390 / 11,372. This single run establishes neither a speedup nor a whole-process bound. Local `wp05-attribution-20261026` artifacts: `native-tests.log`, `phase-regression.log`, `native-tiled256/results.json`, `acceptance.log` and `run_attribution.py`, which reuses the unchanged native runner/provenance and twelve-output comparator. No repeat or benchmark matrix followed.

## Known-size admission increment

On 27 October 2026, `9e24325e29` fixes a shared raw-cache cause of unnecessary scoring rereads rather than adding a calibrated/query payload cache. Previously every miss reserved the full 32 MiB scan maximum, even for a previously decoded immutable scan with a smaller known capacity. The adapter now retains one checked capacity bound per source/visible scan, shared by clones under the existing cache lock. Unknown scans still reserve the maximum. Every decode rechecks its preconversion and retained-capacity requirements, the original hard maximum and oversized-scan diagnostic remain, and a source that exceeds its previous bound fails cleanly. The source must remain immutable. Size metadata uses one `ByteCount` per visible scan in the native metadata/control allowance, not intercepted process-memory accounting. No scoring, calibration, native format, science, scheduling or cap change was made.

| Requirement | Concrete check and observed result |
| --- | --- |
| Avoid unnecessary rereads without a new payload cache | Added one real indexed-file regression in the existing raw-access test. Against unchanged code it fails at the intended unnecessary eviction. After the fix, warm/clear then A/B/A through a clone needs exactly two decodes and one hit, with no eviction, when both measured retained scans fit but a full-maximum reservation would not. |
| Preserve ownership, identity and admission | Exact arrays and acquisition metadata match the native source. Same visible IDs on independent sources do not inherit bounds. Clone sharing, scoring attribution and all ten phase totals reconcile. The same byte cap holds, cancellation rejects new admission, live leases remain readable and last-owner release is zero. |
| Preserve existing failure/public consumer behavior | All five actual vendor-enabled Release tests pass: raw access, raw failure, extraction, scoring and quadratic transformation. The initial test failure exposed a changed oversized-scan diagnostic; that diagnostic was restored without weakening the existing test. |
| Exercise the complete native public workflow | Exactly one unchanged full-library native Bruker tiled256 CLI run exits zero. Both native readers ON, same input/library/scientific settings/threads/private cap and unchanged before/after provenance. All twelve outputs match the saved native normal reference, including finite nonzero IM, IDs/evidence/scores/quantities and the native row counts above. All ten emitted phase fields reconcile, unclassified zero. |

Observed full-run decodes **11,390 to 11,170 (1.93% lower)** and scoring decodes **7,978 to 7,744 (2.93% lower)** relative to the previous single attribution observation. Scoring decoded numeric bytes **57,461,919,792 to 54,711,174,840 (4.79% lower)**. Scoring mutable copies remain **9,432**, copied numeric bytes remain **76,060,458,768**. Parallel interleaving can affect cache counts, so the deterministic native regression establishes the mechanism, not a universal percentage improvement.

Wall **312.19 s versus prior 291.96 s**, peak RSS **1,073,008 versus 1,033,472 KiB**, decoder time **230.754 versus 217.312 s**. These single observations demonstrate **neither a faster tool nor a process-memory improvement**. Private peak 201,325,727 <= 201,326,592 usable bytes, final live zero. The substantial data-access groundwork is not a completed scoring-engine rewrite. Full package and combined acceptance remain open.

Local `wp05-known-bounds-20261027` artifacts: `native-tests-final.log`, retained `native-tests-initial.log`, `native-tiled256/results.json`, `acceptance.log`, `run_known_bounds.py` and `ACCEPTANCE.md`. The pre-fix failure observation is retained in the implementation worker transcript. The run record retains executed revision `3d3ff6b900`; its unpushed commit message alone was corrected to the repository's valid tag, producing published `9e24325e29` with an identical source tree. No repeated workflow or benchmark matrix followed.

## Native-region contract and calibration hotspot

On 27 October 2026, `6495579dde` implements an OpenDIA-private, immutable native column region and allocation-free batch-query contract. It replaces neither the CLI nor the native frame decoder yet. Exact native double m/z, float intensity/point IM, individual RT, identity and separate scalar-mobility/FAIMS annotations are preserved. Capacities are pre-admitted; scan/range copies retain their lease; malformed batches fail before caller output changes. See the [contract, prototype comparison and requirement-level native checks](../native-region-contract.md). Direct Bruker 25-map/75-acquisition and Thermo 301-map/903-acquisition checks pass through the actual source and extractor. The prototype's columnar/batched ideas inform the design, not its fixed-point m/z, mobility collapse or run-wide allocation.

A frozen full native Bruker profile identifies **75.61% of self user-CPU samples** in `SpectrumAccessQuadMZTransforming::getSpectrumById`. Repeated `getMZArray()` calls create/destroy shared owners inside every calibration expression. `c64e647997` obtains the owned array once and iterates values with the identical quadratic/ppm arithmetic. This small shared correction is separate from the private contract and does not bypass calibration, expose mutable cache aliases or change the scorer. Six affected real Release tests pass in 3.94 s, including direct arithmetic/source-isolation/clone/empty-scan regressions.

The unchanged frozen control completes in 309.34 s, user CPU 1,124.07 s, RSS 1,078,984 KiB, all twelve native outputs matching. The first profile has concurrent-library provenance failure. The frozen profile's original observer failed only during expensive DWARF report generation after successful science. A separate recovery verifies unchanged hashes, the same twelve outputs, private cap/final release, zero lost samples and the flat profile without rerunning the engine. Both observer failures remain retained. These are instrumented single observations, not a balanced benchmark or exclusive stage wall times.

**One optimized full CLI observation completes in 231.23 s versus 309.34 s frozen control (25.25% lower wall time); user CPU falls from 1,124.07 to 235.10 s (79.08% lower).** The calibration accessor's self user-CPU share falls from 75.61% to 2.57%. The next largest self symbol is native `stableOrderByKey` at 25.62%, not classifier inference. These are two instrumented observations, not balanced portable benchmark acceptance. RSS increases from 1,078,984 to 1,193,004 KiB, so no memory improvement is established.

The actual complete 4,000-precursor / 24,000-transition native Bruker tiled256 command, science and eight-thread settings are unchanged except output paths. Engine/wrapper exit zero. All twelve outputs match the native normal reference: export 4,699; FEATURE / MS1 / MS2 4,716 each; PRECURSOR 14,409; TRANSITION 28,296; SCORE_MS2 / PEPTIDE / PROTEIN 4,699 / 3,287 / 1,871; precursor / peptide / protein matrices 942 / 942 / 251. IDs, evidence, finite nonzero IM, scores and quantities match under the existing normalization/tolerance rules. All ten phase fields reconcile; unclassified zero. Cache hits/misses/evictions 2,469 / 10,808 / 10,790; private peak 201,326,386 <= 201,326,592 bytes; final live zero. Executable SHA remains unchanged because the correction is in the rebuilt shared core library; actual loaded library paths and before/after library hashes were checked. Local `comparison.json` separates this source revision from the preserved control executable revision.

This was a real observed runtime improvement for the compatibility path, not a speed claim for the then-unintegrated region contract or a quality-matched competitor comparison. The calibrated extraction integration below is now complete. Its first pair does not improve runtime. Cold native index costs, physical-frame reuse, scoring copies, whole-process memory, clean-machine Thermo packaging, full Thermo science and competitive quality remain open. Do not infer that a neural classifier is the present bottleneck from raw-work phase ownership.

## Calibrated region extraction integration

`afba309299` and diagnostic fix `f9f42ff6c6` connect private `OpenDIA::RegionWorkflow` to actual tiled MS1/MS2 extraction. Calibrated native columns retain immutable byte ownership without detached legacy arrays. A generic protected extraction hook and calibration snapshot getters preserve shared users. Other read modes, unsupported wrappers and non-tophat filters keep the reference extractor. Scoring is unchanged and retains wider mobility evidence.

Final acceptance passes 192 public/CLI/installed tests, including 43 installed-consumer steps, twelve verbose boundary tests and six fresh installed consumers. Direct calibrated Bruker/Thermo samples match exact reference chromatograms. Invalid requests, cache calibration identity, native-source isolation, open m/z/IM endpoints, closed RT boundaries, floating-point addition order, admission, cancellation and last-owner release pass. The initial oversize diagnostic failure remains retained; its existing public test passes unchanged after the fix.

One frozen-control/current-candidate pair runs the same complete native Bruker input, 4,000 precursors, 24,000 transitions, eight threads and private cap. Both commands exit zero. All twelve native outputs match. Row counts are export 4,699, FEATURE/MS1/MS2 4,716 each, PRECURSOR 14,409, TRANSITION 28,296, SCORE_MS2/PEPTIDE/PROTEIN 4,699/3,287/1,871 and matrices 942/942/251. Nonzero IM, scores and quantities remain preserved under existing identity/tolerance rules.

| Observed quantity | Frozen control | Region consumer |
| --- | ---: | ---: |
| MS1/MS2 extraction detached copies | 744 / 703 | 0 / 0 |
| Extraction copied numeric bytes | 7,817,201,784 | 0 |
| Scoring detached copies | 9,432 | 9,432 |
| Total decodes | 10,743 | 11,018 |
| Wall seconds | 211.22 | 220.56 |
| User CPU seconds | 216.74 | 223.65 |
| Peak RSS, KiB | 1,244,904 | 1,082,332 |

All ten phase fields reconcile. Unclassified fields are zero. Candidate retained peak is 201,325,702 <= 201,326,592 bytes, with final live zero. Before/after executable/library/input provenance is unchanged within each run. Both readers stay enabled, with no converted intermediate.

**Decision:** accept the compatible extraction integration, not a faster backend. This sequential, unreplicated pair is 4.42% slower with 13.06% lower observed RSS. Separate extraction/scoring payloads can compete for the shared LRU, and scheduling affects cache counts. Do not claim improved physical random access from zero copies. Next test one bounded Bruker physical-frame producer across actual lane requests. Require fewer physical decodes and lower cold-inclusive cost before calling it faster. Full Thermo science and competitor parity remain unverified.

See [requirement-level checks and local evidence](../native-region-contract.md#calibrated-production-consumer-27-october-2026). No repeated benchmark matrix or scorer replacement is part of this increment.

## Next increment

**Superseded implementation restriction, 27 October 2026:** the trace below remains evidence that mutable legacy aliases are unsafe, but its accessor-preserving four-step migration is now only a compatibility option. The owner explicitly authorizes replacing whole access/extraction/scoring components in `src/topp/OpenDIA/`. The calibrated extraction slice above is complete. The next implementation is the bounded physical-frame gate in the [performance strategy](../performance-strategy.md), not mandatory migration through `ISpectrumAccess::SpectrumView`. New scientific scorers use independent quality gates rather than exact old-score equality. Keep all historical observations and failures intact.

### Scoring-boundary trace, 27 October 2026

The trace is complete on `9e24325e29`, but **no scoring-copy reduction has been implemented**. The existing const view is sufficient for extraction, not a drop-in replacement for scoring:

- Every `OpenSwathScoring::fetchSpectrumSwath` overload routes through `ISpectrumAccess::getMultipleSpectra`, which returns `SpectrumSequence`, a vector of mutable `SpectrumPtr`. Even a const `Spectrum` exposes mutable array pointers. Returning cached arrays through that API would break the detached-copy contract.
- `SpectrumAccessQuadMZTransforming::getSpectrumById` changes m/z in place. Its inherited `getSpectrumView` deliberately calls that mutable path. Forwarding the underlying raw view would silently bypass quadratic/ppm calibration. A newly allocated calibrated array would need independently retained byte ownership, including after wrapper destruction.
- DIA integration and `IonMobilityScoring::computeIonMobilogram` consume numerical arrays without modifying them, but their present spectrum input types do not enforce that property. MS1 drift contrast also needs the MS2 sequence. Early filtering to the extraction mobility window would lose wider scoring evidence.
- The existing nested precursor-MS1 reuse already avoids a second fetch in its supported case. Aggregate scoring counters do not identify the MS1/MS2 split, so 9,432 copies is not evidence that 9,432 identical reads can be removed.

In that local-substitution investigation, the master-plan stop rule applied to returning raw aliases through the mutable API. No `const_pointer_cast`, mutable alias of cached columns or speculative query cache was added. This was a code-boundary prerequisite for that adapter, not an external blockage or a prohibition on the now-authorized OpenDIA-owned rewrite.

**Earlier compatibility-migration option, not the mandatory next implementation:**

1. Introduce a read-only scoring path using the existing `ISpectrumAccess::SpectrumView` owners. Start with the single-map, `simple` addition path, retaining existing detached mutable fetch/merge behavior for other paths. Preserve the exact nearest-RT selection/order, full mobility arrays and exception-time release. Keep the existing public mutable API unchanged.
2. Carry that read-only input through the actual DIA and IM numerical consumers, including MS1/MS2 contrast, rather than converting it back into a mutable facade. Keep calibrated wrappers on their correct fallback until a separately charged calibrated view is implemented. Do not claim all scoring copies removed while fallback copies remain.
3. Extend the existing native raw-access regression with the real scoring entry point. Before/after checks must map copy counters to the exercised branch, compare every enabled score and wider-IM evidence, and test mutable isolation, quadratic/ppm fallback, clone/cancellation lifetime and final charge zero. A raw-view microbenchmark alone is insufficient.
4. Build the affected real libraries/tool/tests. After a source change, run one unchanged full native Bruker workflow and compare all twelve outputs, phase counters, private peak and final release. Require observed copy reduction for the exercised end-user path before accepting a performance increment. Report wall time honestly, with no speedup claim from a single observation. Keep both native readers ON and do not convert vendor inputs.

**Boundary verification:** the actual raw-access, quadratic-transform, OpenSwath scoring and ion-mobility tests pass on the unchanged implementation. `DIAScoring_test` was initially not run because its executable had not been built. Building that actual target and running it passes, including its IM scoring sections. Local `wp05-scoring-boundary-20261027` retains the initial not-run result, build output and final DIA test result. These are contract checks, not new copy-reduction or speed evidence. No redundant full native workflow or benchmark matrix followed this documentation-only trace. The previous [whole-stack replay](../VALIDATION.md#current-whole-stack-replay-27-october-2026) and twelve-output acceptance remain the evidence for the implemented stack, not acceptance of this future migration.

Thermo runtime discovery/native workspace attribution and full sharing/calibration acceptance remain separate open gates.

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
