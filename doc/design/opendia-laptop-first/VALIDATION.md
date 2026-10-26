# WP-02 / WP-03 requirement observations

**26 October 2026, 19:10 UTC.** Implementation `4024bc8325d4fecc35a300891ec7937d10c35624`, native Release (`-O3 -DNDEBUG`). This maps the [WP-02 contract](work-packages/02-resource-budget.md) and [WP-03 contract](work-packages/03-raw-tiles.md), including unmet requirements. A passing local primitive or generated-fixture check is not a passing whole-engine gate.

## WP-02

| Requirement | Concrete check and observed result | Scope / remaining gate |
|---|---|---|
| One checked byte budget, reservations and release | `OpenDIAResourceBudget_test`: exact/minimum boundaries, `SIZE_MAX` overflow, invalid configuration/requests, move/split/commit, zero/over-commit, replacement, exceptions and tokens outliving the budget all pass. | Private API passes. Compatibility callers can still omit charges. |
| Cancellation and contending admissions | Same test: per-call stop, permanent cancellation and eight-thread contention pass with released reservations. | No fairness/deadline guarantee or runtime deadlock detector. |
| Bounded queues and payload ownership | `OpenDIABoundedQueue_test`: wraparound/order, move-only payloads, failed-push retention, full/empty stop/close/cancel, reentrant discarded destructors and gated contending producers pass. | Item capacity is not payload-byte accounting. |
| Whole-job acquisition before materialization, downstream drain before producers, no blocking admission while holding a job | `OpenDIAResourceInterop_test` reserves fixture workspace/drain before loading, admits and splits complete jobs, then observes full queues plus blocked producer admission. Writer release drains every row within the test deadline. | Protocol passes in the real-loader fixture pipeline. Production-wide decoder/index/model/evidence/writer admission is **open**, not implied by this test. |
| Preserve immutable record ownership and charge actual retained storage | `OpenDIARecords_test` and resource interop pass shared-copy/final-owner checks, vector/string capacities, long native IDs/provenance, one-byte-undercharged rejection for all three record types, malformed/nonfinite/count-capped records and payload-before-charge destruction. | Validation/by-value copies and later mutable growth need separate admission. |
| Real data survives loader, extractor, queues and writer | Resource interop uses actual `MzMLFile`, public `ChromatogramExtractorAlgorithm`, record constructors and reopened TSV rows. Non-IM **513/513**, IM **32,750/32,750**, FAIMS **44,354/44,354** coordinates/values survive; extraction also matches a direct sum. | Fixture-specific bounded TSV writer, not the production result writer or a streaming decoder. |
| Failure/cancellation wakes users, joins and discards safely | Both queues peak at **2**, their configured capacity. Cancelled IM pipeline writes **0** rows; final reserved/committed charges are **0** after fixture storage destruction. Queue reentrant-destruction checks pass outside-lock discard. | Real CLI error paths also pass below. End-to-end scheduler cancellation belongs to WP-09. |
| Report reserved, committed and queue peaks separately from RSS | Fresh interop reports reserved peak **67,178,672**, committed peak **68,652,976**, simultaneous peak **68,961,200**, usable **68,961,200** bytes; queues **2/2**. Native allowance **134,217,728**, declared total **203,178,928** bytes. | Separate GNU-time RSS in complete CLI runs below. Peaks occur at different instants and are not added. Neither declared allowance is an allocator bound. |
| Oversized work fails explicitly, never truncate populations | Budget invalid/oversized admissions and real `TOPP_OpenDIA_raw_oversize_failure` pass. Full CLI searches retain all **24,000** transitions. Current 128 MiB search fails rather than silently shortening the population. | Splitting/spilling not implemented. Failure is allowed locally, not a successful low-memory search. |
| Measured native reserve and complete live-buffer accounting | Heap profile identifies native metadata/CV/Xerces allocations outside private enforcement. Actual CLI private cache peaks obey caps, but native allocations may precede scan checks. | **Open:** conservative decoder/index/model/downstream bounds, mutable-copy growth, full-process/child/cgroup peaks, final 4 GiB/2 GiB acceptance. A fixed 64 MiB raw allowance does not satisfy this. |
| Disk quotas, platforms and instrumentation | Disk admission is explicitly deferred to WP-08. Linux native/SDK tests pass. Historical resource ASan/UBSan checks pass with the private archive instrumented. | No Windows/macOS, ThreadSanitizer or whole-stack sanitizer claim. Disk/child acceptance has not been run. |

## WP-03

| Requirement | Concrete check and observed result | Scope / remaining gate |
|---|---|---|
| Indexed reader reuse, one scan/tile, factored metadata/numeric columns, no new disk cache | `OpenDIARawAccess_test` exercises the actual indexed reader and returned private records. Exact metadata/arrays match native reads for **114/150/19** scans and **513/32,750/44,354** points, non-IM/IM/FAIMS respectively. | Indexed mzML spectrum maps only. Split-window input unsupported. |
| Separate MS1/MS2; exact RT, m/z, units, frame/scan metadata, isolation and interval edges | RawAccess compares native metadata/arrays and tests long isolation overlap, open/closed RT/isolation edges. Indexed DOM decoder and drift-validation regressions pass standard/nonstandard IM arrays, empty arrays and closed edges. Full converted-input audit preserves 734 spectra, 58,884,702 points and 24 windows. | Converted mzML parity does not validate native per-frame/scan diaPASEF geometry. That complete-engine gate is **open**. |
| Missingness and discrete FAIMS | RawAccess distinguishes absent IM, measured-zero IM, discrete FAIMS zero and empty scans. `TOPP_OpenDIA_raw_faims_failure` prepares a real indexed three-CV fixture and verifies explicit rejection rather than merged lanes. | Lane-aware FAIMS complete-engine support is **open**. Rejection is not support. |
| Preserve public extraction at two limits and increasing query counts | Actual public extractor matches exact native RT/intensity arrays with nonzero references at **1/16/4,096** queries. Fresh Release IM raw peaks are **24,048 / 958,472 bytes**, with **150** misses at every query count for each fixed limit. Output grows **300 to 1,228,800** elements. | Coordinate-grid check only. Outputs, hypotheses and native decoder workspaces are outside this raw-cache bound. |
| Shared clone cache, charges live until last view, separate charged mutable copies | RawAccess passes eviction, retained array/container owners, mutation isolation, cancellation and final-owner release. Copy-admission regression fails on old code and admits the measured fitting copy after repair without unnecessary eviction. Genuine pinned pressure still fails without overrun. | Mutable copies can grow later without interception. |
| Opt-in `readOptions=tiled`, unchanged default/public compatibility | Fresh `TOPPWRITEINI_OpenDIA`, mass/IM defaults, direct/tiled/small non-IM and IM CLI workflows and four parity tests pass. SQLite/Parquet output identity and v0 score handoff tests pass. Existing OpenSwathWorkflow caller/output regressions and installed consumers pass. | No default reader change or public allocator interface added. |
| Explicit invalid/oversized/pinned/cancelled failure; no partial tile | RawAccess checks bad indices/queries, missing index, removed indexed payload and oversized scans/copies. `OpenDIARawFailure_test` returns all **12** indexed-reader admission errors through filtering, both calibrations, legacy and parallel-wave extraction. Real CLI oversize case returns an error, not process termination. | Failure correctness passes. Current 128 MiB complete search still fails extraction and is **not accepted**. |
| Preserve full scientific outputs and IM evidence | Real native OpenDIA direct/tiled workflows preserve four non-IM and six PASEF fixture rows at 65/256 MiB. Full-library TimSim Release searches pass every output check in the next section. Removed-IM negative output check rejects; equal nonfinite values are not accepted. | Tiny-fixture training is too small for scientific score calibration. Complete-library parity closes one converted-input correctness gate only. |
| Reduce memory without assuming speedup | Nine fresh balanced full-library Release processes complete with parity. Tiled256 median RSS decreases **46.45%**, elapsed time increases **2.305x**, below. Larger cache and tiny-fixture counterexamples remain recorded. | No speedup or universal RSS claim. Repeated decoding/shared extraction and whole-job admission remain **open**. |
| Preserve frozen inputs and separate conversion cost | Before/after source/binary/dependency/flags/input/library/comparator/archive hashes match. Complete 2,000-target/2,000-decoy library retained. Conversion separately measured **16.33 s / 299,904 KiB RSS**. | One converted TimSim control run only. Astral and remaining cohort are untested, not declared externally blocked. |
| Native vendor inputs and downstream integration | Current build has Thermo/Bruker readers disabled. Existing separate vendor-enabled FileConverter produced indexed mzML without changing frozen data. | Concrete native-input acceptance constraint. Tiled adapter itself requires indexed mzML. WP-05/06 private scan-view integration and WP-09/10 acceptance remain **open**, not replaced by conversion. |

## Changed native scoring and calibration outputs

`OpenDIARawFailure_test` calls the actual native APIs with a generated indexed pressure fixture, compared with unrestricted native references. This is a precise lifetime regression, not a substitute for the complete CLI searches. Old-core final-charge assertions fail on all nine scoring paths with **1,304 to 2,896** committed bytes. Fresh fixed-core observations report `native_reference=equal committed=0 reserved=0` on every path:

| Changed path | Output / ownership check |
|---|---|
| `calculateDIAScores`: ALL | Fragment and precursor vector plus mass errors are finite and equal within absolute `1e-10`; no retained charge. |
| `calculateDIAScores`: RANKING | Same returned vector comparison for the ranking call; no retained charge. |
| `calculateDIAScores`: REMAINDER | Same returned vector comparison for the remainder call; no retained charge. Existing `OpenSwathScoring_test` separately checks RANKING + REMAINDER against ALL. |
| `calculatePrecursorDIAScores` | Precursor ppm/isotope outputs match; no retained charge. |
| `calculateDIAIdScores` | Identification fragment outputs match; no retained charge. |
| `fetchSpectrumSwath`: access overload | Returned m/z/intensity arrays equal native with positive signal; no hidden scratch owner after returned owners release. |
| `fetchSpectrumSwath`: map overload | Same returned-array and final-charge check. Existing native scoring tests also cover spectrum addition/resampling and IM filtering. |
| Nested MS1 IM reuse | Positive MS1/MS2 drift, weighted drift, deltas, contrast coelution/shape and summed contrast match. Nested owners survive until outer contrast finishes, then release. |
| Exception exit | Actual native `ParseError` survives after both reads; final charge is zero. |
| `correctIM` and `correctMZ` | Both actual calibration APIs preserve fit inputs/results and positive windows within `1e-10`; usable = peak **2,512 bytes**, final charge **0**. Both pressure cases fail on old core. |

The scoring vector includes mass deviation, weighted mass deviation, isotope correlation/overlap, DIA dot product/Manhattan, b/y series, precursor ppm/isotope scores and mass errors. No claim that every branch of every caller has a dedicated pressure fixture is made. The real caller regressions and complete workflow provide separate integration coverage.

## Complete end-user outputs and resource observations

The project's `src/topp/OpenDIA/test/compare_results.py` reopened actual complete CLI results. It checks identities exactly, numeric values at relative `1e-6` / absolute `1e-8`, finite nonzero IM and all feature/score outputs. Only generated run/file/feature identifiers are normalized. Fresh replay also checked all three matrices:

| Public output | Observed rows / columns, equal across complete workflows |
|---|---:|
| `results.tsv` | 4,688 / 55 |
| `FEATURE` | 4,700 / 10 |
| `FEATURE_MS1` | 4,700 / 18 |
| `FEATURE_MS2` | 4,700 / 38 |
| `FEATURE_PRECURSOR` | 14,219 / 4 |
| `FEATURE_TRANSITION` | 28,200 / 43 |
| `SCORE_MS2` | 4,688 / 6 |
| `SCORE_PEPTIDE` | 3,287 / 6 |
| `SCORE_PROTEIN` | 1,871 / 6 |
| Precursor matrix | 985 / 6 |
| Peptide matrix | 985 / 3 |
| Protein matrix | 251 / 2 |

Nine balanced complete processes, three per reader, use unchanged full input/library/science, eight threads, warm filesystem cache and a shared host. The [raw sample table](work-packages/03-raw-tiles.md#release-completed-scoring-scans) remains available. Median normal: **162.44 s / 1,834,100 KiB RSS**. Tiled256: **374.40 s / 982,188 KiB**. Tiled1024: **173.53 s / 1,852,132 KiB**, no clear benefit. All private peaks stay within **201,326,592 / 1,006,632,960 bytes**, respectively, with final live-copy bytes zero. This is a measured memory/time tradeoff, not global budget enforcement. The earlier failed matrix remains rejected.

The original failed256 command now completes **349.94 s / 1,104,268 KiB**, with all outputs matching. That focused success is not a timing comparison. Current128 on the same published Release source completes calibration but fails extraction admission: engine and wrapper **exit 8**, **78.74 s**, no final exports. Its partial RSS is excluded from completed-run comparisons. The former Debug128 failures remain historical evidence, not a controlled repair comparison. The smaller cache also failed to save process RSS in the separate tiny-fixture matrix.

## Replay and acceptance decision

Fresh actual project replay at **19:08-19:10 UTC** passed **204/204** selected native class/CLI/SDK tests, **6/6** nested installed consumers, **7/7** verbose public-interface tests and the full real-output replay. Per-requirement observations above, not those aggregate counts, support each claim. Earlier standalone Debug/NDEBUG pass **10 each**; private raw ASan/UBSan pass **2**, with leak detection off and the native SDK not instrumented.

Reproduce the focused public-interface checks in the actual native build, plus the real CLI output/consumer fixtures:

```bash
ctest --test-dir ../OpenMS-wp03-release-build -R '^(OpenDIAResourceBudget_test|OpenDIABoundedQueue_test|OpenDIAResourceInterop_test|OpenDIARawAccess_test|OpenDIARawFailure_test|OpenSwathScoring_test|SwathMapMassCorrection_test)$' -V
ctest --test-dir ../OpenMS-wp03-release-build -R 'OpenDIA|OpenSwathWorkflow|TestExternalCode' --output-on-failure -j3
ctest --test-dir ../OpenMS-wp03-release-build/installed-consumer/build --output-on-failure
```

Local command vectors, test selection, native/nested/verbose logs and reopened-output observations are retained in `representative-1312/acceptance-recheck-1908`; full process evidence remains in `scoring-balanced-1740`, `scoring-acceptance-1740` and `current-128-1900`. Raw evidence sanitization/publication is still a follow-up, so these local artifacts are not presented as downloadable public evidence.

**Decision:** the local lifetime repair and one-input full-output parity are verified. WP-02/WP-03 as complete packages are **not accepted**. Retain the native-geometry, broader-cohort, repeated-decode, 128 MiB, global-accounting and downstream acceptance gates above. No source change, cap inflation, scientific truncation, merge or issue closure is justified by this replay. Keep branches stacked and obtain owner approval only after combined-design validation.
