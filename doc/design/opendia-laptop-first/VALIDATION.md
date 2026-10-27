# WP-02 / WP-03 requirement observations

**Historical snapshot, 26 October 2026.** Implementation `4024bc8325d4fecc35a300891ec7937d10c35624`, native Release (`-O3 -DNDEBUG`). The tables below map the [WP-02 contract](work-packages/02-resource-budget.md) and [WP-03 contract](work-packages/03-raw-tiles.md) at that revision, including unmet requirements. They are not claims that later vendor/read-only work is absent. Current native/read-only and known-size admission results are in [WP-05](work-packages/05-shared-extraction.md); the prior whole-stack replay is recorded below. A passing local primitive or generated-fixture check is not a passing whole-engine gate.

## Latest scoped validation, 27 October 2026

After the whole-stack replay below, `6495579dde` adds the private native column contract and `c64e647997` fixes shared calibration array-owner churn. The actual vendor-enabled Release tool and six affected tests rebuild and pass. Direct region-to-source/extractor comparisons pass for Bruker 25 maps / 75 sampled acquisitions and Thermo 301 maps / 903 sampled acquisitions, with exact identity/columns, edge/failure/lease checks mapped in the [contract observations](native-region-contract.md#observed-requirement-to-check-mapping-27-october-2026).

The full optimized native Bruker CLI preserves all twelve scientific outputs, reconciles every phase counter, respects its private cap and ends with zero live charge. Observed wall 309.34 to 231.23 s, user CPU 1,124.07 to 235.10 s, RSS 1,078,984 to 1,193,004 KiB. See [full results and retained observer failures](work-packages/05-shared-extraction.md#native-region-contract-and-calibration-hotspot). These are instrumented single observations, not a balanced speed or lower-RSS guarantee. The private region is not yet consumed by the CLI.

**Final-revision replay:** the complete original 188-test public/native/installed selection was rerun on `c64e647997`, together with the new region and DIA/IM tests: **191/191 pass**, including **43 actual installation/consumer steps**. An additional **11/11 verbose boundary replay**, **6/6 freshly installed consumers** and both direct native region paths pass. All twelve saved full native scientific outputs were reopened and compared again, without launching another performance workflow. New contract headers remain private, with no new installed API or dependency.

| Requirement / changed output | Final-revision replay observation |
| --- | --- |
| WP-02 admission, ownership and pipeline | Budget/queue overflow, stop, cancellation, exception and contention tests pass. The real loader/queue/writer reopens 513 / 32,750 / 44,354 rows exactly; cancelled IM writes zero. Both queue peaks are 2, simultaneous accounted peak 68,961,200 bytes equals the usable cap, final-release assertions pass. |
| WP-03 exact source/extraction behavior | Metadata and arrays match for 114 / 150 / 19 scans and 513 / 32,750 / 44,354 points. Actual extraction at 1 / 16 / 4,096 queries passes at both caps, with and without point IM where supported. |
| Calibration and scoring integration | Both IM/non-IM calibration comparisons report native reference equal and final charge zero. All nine scoring ownership paths report equal reference and committed/reserved zero. Quadratic/ppm arithmetic, clone/source isolation, empty scans and unchanged intensity tests pass. |
| New bounded region contract | Endpoint/order/duplicate, missingness/unit, capacity/overflow, cancellation, immutable-view lifetime and invalid-batch untouched-output assertions pass on final source. Direct Bruker 25 maps / 75 acquisitions and Thermo 301 maps / 903 acquisitions pass through the native source, region and actual extractor: 326 three-acquisition regions checked. |
| Public CLI and failures | All selected direct/tiled IM/non-IM workflows and parity comparisons pass, including SQLite/Parquet identity, INI/CTD/defaults, explicit oversize/FAIMS failures and OpenSwathWorkflow regressions. Reopened full Bruker export, eight feature/score tables and three matrices match all identities/evidence/scores/quantities and finite nonzero IM under the established normalization/tolerance rules. |
| Local packaging boundaries | Full/core/relocated/no-Boost install/configure/build/run steps pass. Fresh installed Qt, CLI, OSWParquet, OpenDIA, core and PeptDeep consumers pass. This is a local Linux package replay, not cold clean-machine vendor-runtime distribution acceptance. |
| Measured improvement and limits | Frozen/optimized library hashes and phase/cap/release observations recheck successfully. The 25.25% lower observed wall and 79.08% lower user CPU belong to the shared calibration correction, not the unintegrated new region. Higher RSS remains recorded. No additional benchmark or broader quality claim. |

Evidence: `wp05-data-contract-20261027/final-stack-replay/{selection.json,build.log,tests.log,boundaries.log,installed.log,native-regions.log,output-recheck.log,profile-provenance-recheck.log,requirement-observations.json}` and `replay_final.py` / `check_final_replay.py`. Default calibrated region/scoring integration, cold index/workspace cost, physical-frame reuse, automatic 128 MiB admission, process-memory bounds, full Thermo science/runtime packaging and competitor-quality acceptance remain open.

<a id="current-whole-stack-replay-27-october-2026"></a>

## Prior whole-stack replay, 27 October 2026

Implementation `9e24325e29`, including all earlier stacked WP-01/02/03/05 changes, was rebuilt before this replay. Actual project native/CLI/installed-package checks pass **188/188**, including **43 real installation/consumer steps** and the 800-header no-Boost consumer build. A separate verbose boundary replay passes **7/7**, and the freshly installed consumers pass **6/6**. The requirement observations, rather than these totals alone, are:

| Requirement / changed public output | Current executed check and observed result |
| --- | --- |
| WP-02 checked admission, overflow and cancellation | Native resource-budget tests pass minimum/exact/oversized/overflow, move/split/commit, replacement/exception, stopped/cancelled and contending admissions. Queue tests pass ordering, full/empty stop/close/cancel, failed-push ownership and reentrant destruction. |
| Whole real-loader/extractor/queue/writer path | Resource interop reopens **513/513 non-IM, 32,750/32,750 IM and 44,354/44,354 FAIMS rows**, exactly. Both queues peak at **2**; simultaneous accounted peak **68,961,200 <= 68,961,200 usable bytes**. Cancelled IM writes **0** rows and the test requires final reservations/charges zero. This remains a fixture pipeline, not production-wide workspace admission. |
| WP-03 metadata, numeric columns, units/missingness and extraction axes | Actual raw-access test matches native metadata and arrays for **114/150/19 scans**, **513/32,750/44,354 points**, and public extraction at **1/16/4,096 queries** at both cache limits. Its interval-edge, missing/zero mobility, FAIMS-zero, empty-scan and invalid-query checks pass. |
| Cache/copy/calibrated ownership, including current change | Native clone, pinned/admission, mutation-isolation, cancellation and independent array/container final-owner checks pass. Read-only IM/non-IM extraction and quadratic/ppm fallback pass. Known-size A/B/A observes **two decodes plus one hit**, no unnecessary eviction, unchanged cap and exact arrays; source isolation and ten-field phase reconciliation pass. |
| Changed scoring/calibration outputs and failures | All **nine** scoring paths report `native_reference=equal committed=0 reserved=0`: ALL, RANKING, REMAINDER, precursor, ID, both fetch overloads, MS1-IM reuse and exception. IM/non-IM calibration reports equal native reference and final charge zero. All twelve native admission-error paths pass, with the original oversized-scan CLI diagnostic preserved. |
| Public CLI/defaults and output compatibility | Actual direct/tiled/small non-IM and IM workflows plus four parity comparisons pass, as do oversize/FAIMS explicit failures, SQLite/Parquet identity, both v0 score handoffs, INI/CTD/help/default checks and existing OpenSwathWorkflow regressions. FAIMS rejection is not lane-aware support. |
| Real packaging and integration | Current libraries/headers are installed, configured, compiled and run through full/core/relocated consumers. Core-only correctly rejects CLI discovery; relocated consumers run; no-Boost compiles all **800** public-header translation units. Fresh installed Qt, CLI, OSWParquet, OpenDIA, core and PeptDeep consumers all pass. This is Linux acceptance, not clean-machine vendor-runtime distribution acceptance. |
| Complete native scientific result and resources | The same current source tree already completed the full native Bruker CLI and all **twelve** native-reference output comparisons documented in [WP-05](work-packages/05-shared-extraction.md#known-size-admission-increment), including finite IM, unchanged full library/science and stable provenance. Current-tree equality is checked after the message-only commit amendment. Private peak **201,325,727 <= 201,326,592**, final live zero, all phase fields reconcile and unclassified zero. No second complete-library run was needed. |
| Layout/performance outcome, not an inferred improvement | The factored production scan path, read-only extraction and detached mutable scoring boundary are exercised above. The prior controlled layout comparison is historical, not rerun or converted into an isolated engine-speed claim. Latest scoring decodes **7,978 to 7,744**, but copies remain **9,432**; wall **291.96 to 312.19 s** and RSS **1,033,472 to 1,073,008 KiB**. Thus **no faster tool, whole-process memory improvement or isolated production AoS/SoA gain is established**. |

Local evidence: `wp05-known-bounds-20261027/whole-stack-{build,tests,observations,installed}.log`, alongside the unchanged full-run `native-tiled256/results.json`, `acceptance.log` and `ACCEPTANCE.md`. The fresh replay observes the whole currently implemented result, not just the last regression. Full Thermo science/runtime packaging, broader acquisition/cohort coverage, automatic 128 MiB scheduling, whole-job/native workspace/disk/child accounting and final speed/memory acceptance remain **open and untested here**, not externally blocked or accepted. WP-02/WP-03 remain in review, WP-05 in progress, with no merge or package closure.

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
| Oversized work fails explicitly, never truncate populations | Budget invalid/oversized admissions and real `TOPP_OpenDIA_raw_oversize_failure` pass. Full CLI searches retain all **24,000** transitions. Current automatic-scheduling 128 MiB search fails rather than silently shortening the population. | Splitting/spilling not implemented. Failure is allowed locally, not a successful low-memory search. |
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
| Explicit invalid/oversized/pinned/cancelled failure; no partial tile | RawAccess checks bad indices/queries, missing index, removed indexed payload and oversized scans/copies. `OpenDIARawFailure_test` returns all **12** indexed-reader admission errors through filtering, both calibrations, legacy and parallel-wave extraction. Real CLI oversize case returns an error, not process termination. | Failure correctness passes. Current automatic-scheduling 128 MiB complete search still fails extraction and is **not accepted**. |
| Preserve full scientific outputs and IM evidence | Real native OpenDIA direct/tiled workflows preserve four non-IM and six PASEF fixture rows at 65/256 MiB. Full-library TimSim Release searches pass every output check in the next section. Removed-IM negative output check rejects; equal nonfinite values are not accepted. | Tiny-fixture training is too small for scientific score calibration. Complete-library parity closes one converted-input correctness gate only. |
| Reduce memory without assuming speedup | Nine fresh balanced full-library Release processes complete with parity. Tiled256 median RSS decreases **46.45%**, elapsed time increases **2.305x**, below. Larger cache and tiny-fixture counterexamples remain recorded. | No speedup or universal RSS claim. Repeated decoding/shared extraction and whole-job admission remain **open**. |
| Preserve frozen inputs and separate conversion cost | Before/after source/binary/dependency/flags/input/library/comparator/archive hashes match. Complete 2,000-target/2,000-decoy library retained. Conversion separately measured **16.33 s / 299,904 KiB RSS**. | One converted TimSim control run only. Astral and remaining cohort are untested, not declared externally blocked. |
| Native vendor inputs and downstream integration | At the historical WP-03 revision both vendor readers were disabled and the adapter required indexed mzML. **Superseded by WP-05:** both readers are now enabled, original Bruker `.d` / `.d.zip` and Thermo `.raw` are read directly, native-reader checks pass, and the full native Bruker normal/tiled twelve-output comparison passes. No indexed-mzML conversion is required or performed. | Full Thermo science, clean-machine runtime packaging, native workspace accounting, lane-aware FAIMS and broader WP-09/10 acceptance remain open. See [current evidence and retained cross-representation failure](work-packages/05-shared-extraction.md#native-vendor-input-increment). |

## Layout decision through the production workflow

The [production scan-view follow-through](experiments/layout/RESULTS.md#production-scan-view-workflow-follow-through) maps metadata factoring, numeric IM columns, mutable compatibility ownership, complete scientific outputs and process costs to the actual adapter/public interfaces above. It connects the original layout experiment to nine completed native Release searches without presenting reader/cache/copy/scheduling effects as an isolated AoS/SoA improvement. The v0 `RawTile` contract is unchanged and the legacy bridge still deep-copies numeric arrays.

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

## Low-cap scheduling hypothesis

The [128 MiB counterfactual](work-packages/03-raw-tiles.md#scheduling-counterfactual-at-128-mib) exercises the actual complete native CLI, not a substitute pipeline. The frozen original automatic-scheduling command remains failed. Only output/work directories and the existing outer-loop thread knob differ:

| Requirement / hypothesis | Concrete observation | Acceptance limit |
|---|---|---|
| Is concurrency contributing to active scan pressure? | Actual unchanged-command GDB stack reaches MS1 extraction through the mutable m/z-transform adapter while seven peers are in MS1 transform/extraction. With outer-loop threads set to 1, the native engine completes all 24,000 transitions. | Supports live-copy concurrency as one pressure source. No universal scheduler/job-size bound. GDB's exact-text assertion failed on a truncated string, so that diagnostic is not engine completion. |
| Preserve full input, cap and science | Command-vector comparison permits only output/work paths and `outer_loop_threads -1 -> 1`; global eight threads, 128 MiB total, 32 MiB scan maximum and full 2,000-target/2,000-decoy library remain unchanged. Before/after provenance matches the nine-run matrix. | Explicit scheduling differs from the original automatic command. |
| Complete actual downstream outputs | Native engine GNU-time exit 0, completion log and recovered service journal agree. Project comparator reopens every output in the table above, including finite nonzero IM and all three matrices, with full parity. | The wrapper was killed after 600 s and returns 124. This is recovered native scientific completion, **not a clean timing trial** or a successful wrapper exit. |
| Retain charges within this private cap | Peak 67,108,737 <= usable 67,108,864 bytes, final live 0; 1,106 hits / 12,139 misses / 12,134 evictions. Native observed 1,054.84 s / 760,764 KiB RSS. | Do not add this descriptive recovered observation to the balanced performance matrix or infer a 128 MiB process bound. |
| Reuse the smallest working control | The existing `outer_loop_threads=1` completes this measured configuration, with no production source change, raised cap or scientific truncation. | Automatic admission, conservative whole-job workspace bounds and shared extraction remain **open**. `maxConcurrentSwaths` does not limit this tiled/on-disc path. |

Original supervisor evidence stays preserved in `current128-outer1-1923/results.json`, with explicit recovery status in `recovery-results.json`; the original unchanged-command failure remains in `current-128-1900`. The observation does not override either failure or accept WP-02/WP-03 as complete packages.

### Clean explicit-scheduling repeat

The fresh `current128-outer1-clean-2003` run closes the interrupted-supervisor gap for the explicit scheduling configuration, not the original automatic command:

| Requirement | Concrete observation | Acceptance limit |
|---|---|---|
| Complete a real process without recovery | Independent bounded supervisor runs the same native CLI. Only output/work paths differ from the earlier explicit-outer-1 command. Native engine and native-service wrapper both exit **0**; strict original process/GNU-time checks pass. | Not automatic admission. The earlier wrapper **124** remains rejected. |
| Preserve every public output | Independent project comparator reopens all twelve outputs in the table above with exact scientific identities, declared numeric tolerance, finite nonzero IM and all three matrix levels equal to the normal reference. | Same complete input/library/science, not broader-cohort or native-vendor acceptance. |
| Preserve private cap and ownership | Peak **67,108,737 <= 67,108,864 bytes**, final live **0**. Cache **1,098 hits / 12,147 misses / 12,142 evictions**. Before/after provenance equals the nine-run matrix. | Raw cache/copy accounting, not complete live-workspace enforcement. |
| Observe process and service costs separately | GNU time **1,034.77 s / 857,292 KiB RSS**. Exact native-service invocation journal records **849,977,344 bytes** cgroup memory peak and **0** swap; captured live guards retain **1,200 s / 12 GiB / zero swap / 800% CPU**. | One clean process sample, not a balanced performance comparison, measured universal reserve or 128 MiB process bound. Do not add it to the nine-run matrix. |

The first read-only observer failed because systemd had collected the finished service and returned default exit properties. Its failed log is retained as `observer-first-failed.log`. The corrected check does not accept those defaults: it verifies the actual stored wrapper exit, strict complete engine record, captured live guards, exact invocation journals and reopened scientific outputs. The independent supervisor and native run did not fail. The original automatic-scheduling engine/wrapper **8** remains failed with no exports. Fresh **20:22 UTC** checks also pass **7/7** actual native public-interface tests and **6/6** installed consumers, with per-path ownership observations retained in the verbose log.

## Replay and acceptance decision

Fresh actual project replay at **19:08-19:10 UTC** passed **204/204** selected native class/CLI/SDK tests, **6/6** nested installed consumers, **7/7** verbose public-interface tests and the full real-output replay. Per-requirement observations above, not those aggregate counts, support each claim. Earlier standalone Debug/NDEBUG pass **10 each**; private raw ASan/UBSan pass **2**, with leak detection off and the native SDK not instrumented.

Fresh continuation replay at **19:46 UTC** again passes **7/7** actual native public-interface checks and **6/6** installed consumers. It independently reopens the recovered single-outer-job outputs, verifies exact scheduling-only command differences, all twelve output comparisons, private caps/final release and unchanged provenance. The strict original process checker still rejects wrapper exit124, and the original automatic-scheduling exit8 remains failed. The layout workflow map separately matches the nine completed balanced commands, medians, output counts, public extraction observations and scope limits.

Reproduce the focused public-interface checks in the actual native build, plus the real CLI output/consumer fixtures:

```bash
ctest --test-dir ../OpenMS-wp03-release-build -R '^(OpenDIAResourceBudget_test|OpenDIABoundedQueue_test|OpenDIAResourceInterop_test|OpenDIARawAccess_test|OpenDIARawFailure_test|OpenSwathScoring_test|SwathMapMassCorrection_test)$' -V
ctest --test-dir ../OpenMS-wp03-release-build -R 'OpenDIA|OpenSwathWorkflow|TestExternalCode' --output-on-failure -j3
ctest --test-dir ../OpenMS-wp03-release-build/installed-consumer/build --output-on-failure
```

Local command vectors, test selection, native/nested/verbose logs and reopened-output observations are retained in `representative-1312/acceptance-recheck-1908`; full process evidence remains in `scoring-balanced-1740`, `scoring-acceptance-1740` and `current-128-1900`. Raw evidence sanitization/publication is still a follow-up, so these local artifacts are not presented as downloadable public evidence.

**Decision:** the local lifetime repair and one-input full-output parity are verified. WP-02/WP-03 as complete packages are **not accepted**. Retain the native-geometry, broader-cohort, repeated-decode, automatic 128 MiB scheduling, global-accounting and downstream acceptance gates above. No source change, cap inflation, scientific truncation, merge or issue closure is justified by this replay. Keep branches stacked and obtain owner approval only after combined-design validation.
