# Measured raw-point layout tradeoffs

Measured October 25, 2026. Planning-only experiment, not a production migration.
Protocol, code and reproduction commands: [README](README.md).
All 35 cases, 105 processes, individual timing samples and hashes: [results](results-20261025.json).

## Decision

**Factor acquisition metadata first. Do not declare SoA a universal winner.**

- Factoring metadata out of each point reduced retained layout allocation by **70-92%** in the tested large cases. Flat SoA alone saved only about **2-3%**.
- For numeric streaming, the factored SoA was **5.85-9.08x faster than direct full-record AoS** at the large scales, excluding by-value accessor overhead. It was **2.43-3.97x faster than factored AoS**.
- That does **not** translate into the same extraction speedup. Existing public extraction includes scan traversal, allocations, optional-IM adaptation and query work. Factored SoA was strongly favorable for FAIMS without continuous IM, but only modestly faster than full-record AoS for the tested IM case, and slower than factored AoS there.
- Sparse numeric gathers and complete-row reads can favor AoS. A numeric scan view should be the normal hot path, with explicit row reconstruction when metadata is actually needed. Callers should not manually join a dozen columns.

A reasonable next production design is acquisition metadata plus offsets and an owning numeric scan/tile view, not an immediate conversion of every record family to SoA. Prefer contiguous numeric columns for no-IM streaming. Keep the IM choice open until the real future extractor adapter is implemented and measured. No production files or v0 contract behavior were changed by this experiment.

## What was compared

| Name | Storage/access |
| --- | --- |
| `rawtile` | Actual unchanged WP-01 RawTile, public by-value point accessor |
| `aos` | Same full RawPoint payload, direct const-reference access for numeric loops |
| `soa` | Separate columns, still repeating all per-point metadata |
| `factored-aos` | Acquisition metadata once, offsets, 40-byte hot numeric point records |
| `factored-soa` | Same acquisition metadata and offsets, separate numeric columns |

Every representation also owns an identical acquisition directory and offsets, including empty scans. Thus the RawTile baseline here is **RawTile plus that directory**, not bare RawTile. Full-row gathers return equivalent complete values for every variant, including copying native-ID strings. The direct AoS baseline isolates numeric accessor-copy overhead, not all possible row-copy overhead.

On this ABI, `sizeof(RawPoint)=136`, `sizeof(Coordinates)=128`, `sizeof(HotPoint)=40` and `sizeof(ScanMetadata)=120`. Native-ID allocations and allocator overhead are additional. Factored SoA omits the whole mobility column if all points lack continuous mobility. When mobility exists, it retains `optional<Mobility>` per point, including unit and presence, so both factored variants use about 40 bytes per point before acquisition overhead. This is not a packed validity bitmap or a scan-shared IM-unit implementation.

## Retained memory

Median live allocation delta after layout construction, in bytes per point. This is glibc allocator accounting, **not RSS or whole-engine memory**. The same canonical rows and loaded mzML are retained for every variant, and are excluded by the before/after delta.

| Fixture and measured rows | RawTile | Direct AoS | Flat SoA | Factored AoS | Factored SoA |
| --- | ---: | ---: | ---: | ---: | ---: |
| Non-IM, 1,049,085 | 164.67 | 164.67 | 160.67 | 68.67 | 44.67 |
| FAIMS, 1,064,496 | 200.08 | 200.08 | 196.11 | 40.08 | 16.09 |
| IM, 1,080,750 | 136.59 | 136.59 | 132.60 | 40.59 | 40.59 |
| IM, 4,224,750 | 136.59 | 136.59 | 132.60 | 40.59 | 40.59 |

At 4,224,750 IM points, retained allocation fell from **550.33 MiB to 163.55 MiB**, a **70.28%** reduction, for either factored layout. At 1,064,496 FAIMS points, factored SoA used **16.33 MiB versus 203.12 MiB**, a **91.96%** reduction. FAIMS is acquisition metadata, not a reason to allocate a continuous IM column.

The non-IM fixture has only 513 points across 114 acquisitions. Repetition produces 233,130 acquisitions for its million-point case, so per-acquisition metadata remains substantial. The IM fixture has 32,750 points across 150 acquisitions, and FAIMS has 44,354 points across 19 acquisitions. These are materially different acquisition shapes, not interchangeable estimates of bytes per point for all DIA data.

The JSON also retains baseline RSS, post-construction RSS/peak and final process peak. Validation subsequently creates an independent full AoS oracle, so final peak must not be presented as the candidate layout footprint. Build time is diagnostic only: actual RawTile performs production validation, experimental packers do not repeat the same full validation.

## Speed

Milliseconds per workload, median of three independent process medians. Each process warms up then records seven samples. Lower is better. Dense reads every point. Sparse numeric/full-row workloads visit up to 4,096 deterministically selected row IDs. Extraction uses 34 sorted targets per acquisition, with 0.05 Th narrow and 10 Th wide windows and an optional 0.1-wide IM filter.

### Dense numeric streaming

| Fixture and rows | RawTile | Direct AoS | Flat SoA | Factored AoS | Factored SoA |
| --- | ---: | ---: | ---: | ---: | ---: |
| Non-IM, 1,049,085 | 15.799 | 7.475 | 1.131 | 2.876 | 1.156 |
| FAIMS, 1,064,496 | 54.482 | 10.834 | 1.147 | 4.735 | 1.193 |
| IM, 4,224,750 | 63.841 | 30.488 | 4.905 | 12.692 | 5.215 |

The RawTile-to-direct-AoS difference is primarily an accessor/copy-path difference, not evidence that AoS itself is that much slower. The direct-AoS-to-SoA comparison removes that confounder for numeric loops.

### Actual existing public extraction kernel, including layout adaptation

AoS variants pack m/z/intensity vectors per scan. SoA variants pass their actual numeric vectors directly. Optional IM still needs adaptation in all variants. Cursors advance across sorted targets as in the existing extractor. These are actual `extract_value_tophat` calls, not an experimental substitute.

| Fixture | Layout | Narrow ms | Wide ms |
| --- | --- | ---: | ---: |
| Non-IM, 1,049,085 points | RawTile | 104.570 | 133.479 |
| | Direct AoS | 94.564 | 124.740 |
| | Flat SoA | 76.073 | 103.842 |
| | Factored AoS | 87.066 | 115.417 |
| | Factored SoA | 82.188 | 110.336 |
| FAIMS, 1,064,496 points | RawTile | 64.996 | 83.155 |
| | Direct AoS | 14.176 | 23.074 |
| | Flat SoA | 1.415 | 6.860 |
| | Factored AoS | 11.078 | 18.301 |
| | Factored SoA | 1.371 | 6.772 |
| IM, 4,224,750 points | RawTile | 174.654 | 374.239 |
| | Direct AoS | 131.669 | 330.904 |
| | Flat SoA | 95.238 | 296.967 |
| | Factored AoS | 115.945 | 315.173 |
| | Factored SoA | 128.622 | 326.426 |

At the large IM scale, factored SoA is only **1.024x/1.014x** faster than direct full-record AoS for narrow/wide public extraction, despite its **5.85x** dense-streaming advantage. Factored AoS is faster than factored SoA for those extraction paths while retaining effectively the same memory. Flat SoA has the best observed IM public-kernel medians here, but retains **534.23 MiB** rather than about **163.55 MiB**. No single candidate wins speed, memory and all access patterns.

### Sparse/full-row counterexamples and uncertainty

- Large IM sparse numeric: direct AoS **0.0446 ms**, factored AoS **0.0369 ms**, factored SoA **0.0641 ms**. Reading related values from one compact row can beat separate columns.
- Large IM full-row reads: direct AoS **0.3955 ms**, factored AoS **1.0253 ms**, factored SoA **1.2149 ms**. Factored gathers perform an offsets lookup and reconstruct metadata. A scan-aware caller could avoid repeated lookup, but that alternative was not measured.
- Large non-IM full-row reads: direct AoS **0.2879 ms** versus factored SoA **1.3688 ms**, about **4.75x slower**. This cost must not disappear from the design decision.
- Large IM factored-AoS public narrow process medians ranged **115.80-118.11 ms**, versus factored-SoA **125.94-174.49 ms**. Dense ranged **12.66-13.06 ms** versus **4.90-7.00 ms**, respectively. Both the streaming win and the extraction counterexample are visible across the three trials.
- The host was not exclusive. Some trials were noisy: FAIMS RawTile public narrow ranged **37.32-160.76 ms**, and non-IM factored-SoA public wide ranged **109.98-310.64 ms**. All samples are retained, no outliers were dropped. Treat precise ratios and small differences as provisional, especially when ranges overlap. There are no confidence intervals or statistically established universal rankings.

Separate experimental lower-bound extraction results are retained in JSON as `narrow` and `wide`. They are supporting access-pattern evidence only. In particular, the large IM factored-SoA experimental wide loop was slower than factored AoS, another reason not to infer extraction performance from dense bandwidth alone.

## Reproduction and acceptance observations

- Production source: `3cd24a25026615b26671a46596b3eeae7de89986`, unchanged and clean.
- Final measurement code: planning commit `81caeb391b`, with complete source/binary/fixture SHA-256 values in JSON.
- Measurement interval: **18:41:54-18:50:13 UTC, October 25, 2026**.
- Linux x86-64, Intel Core i7-1370P, CPU **10** pinned, sibling **10-11**, **24 MiB** reported LLC, `powersave` governor unchanged. Other host activity was not controlled.
- GCC **13.3.0**, C++23, `-O3 -DNDEBUG` for all layouts and the unchanged production records and extractor translation units. Both actual public extractor overloads and `extractChromatograms` were verified as executable-defined symbols. The remaining shared OpenMS SDK is Debug with STL debug mode off, so these are **not whole-engine Release results**.
- **105/105** fresh processes succeeded: five variants, three fixtures at native and million-point repeated sizes, plus four-million-point IM, three randomized trials, seed `20261025`, seven timed repetitions after warm-up.
- All full-field row roundtrips, acquisition offsets and empty scans were checked. Missing versus zero mobility/FAIMS, units, RT, isolation, IDs and native strings were preserved. Real-file extraction matched the actual public `extractChromatograms` interface for m/z and supported per-point IM subsets. Every timed public-kernel result was also compared individually with the reference, not only by aggregate checksum.
- **21/21 Release CTests** passed, including all five layouts on four real mzML fixtures and always-active edge checks for empty/singleton scans, duplicate m/z, exact boundaries, mixed/missing mobility units, explicit FAIMS zero and long native IDs.
- Runner self-checks passed. The result validator accepted the full real matrix, matched all cross-layout/trial workload checksums, and independently rejected 12 corrupted/incomplete result cases. These reporting tests are not a substitute for public-interface correctness.
- Existing native/public CLI and installed-SDK regression paths also passed **50/50** and **6/6**, respectively. This guards the unchanged implementation and integration environment, not production acceptance of an unimplemented layout migration.

## Requirement-to-observation closure

A follow-up on **October 25, 2026 at 18:54-18:55 UTC** exercised the native project path as well as rechecking the experimental evidence. The actual native targets `ChromatogramExtractorAlgorithm_test`, `ChromatogramExtractor_test`, `SpectrumAccessOpenMS_test`, `MzMLFile_test` and `OpenSwathChromatogramExtractor` were built successfully. The tests and CLI use the actual shared `libOpenMS`, verified by linkage, not the benchmark's optimized extractor translation unit.

```bash
cmake --build <native-build> --target ChromatogramExtractorAlgorithm_test \
  ChromatogramExtractor_test SpectrumAccessOpenMS_test MzMLFile_test \
  OpenSwathChromatogramExtractor -j3
ctest --test-dir <native-build> \
  -R '^(ChromatogramExtractorAlgorithm|ChromatogramExtractor|SpectrumAccessOpenMS|MzMLFile)_test$|^TOPP_OpenSwathChromatogramExtractor' \
  -V -j2
```

**Observed: 14/14 passed in 6.48 seconds.** This includes four native class tests and five genuine extraction CLI runs, each with its existing FuzzyDiff output comparison. Every generated chromatogram file matched its expected project fixture. Nonblocking existing INI version warnings and an existing range-loop compiler warning were retained in the logs, not treated as failures. These are native regression observations, not timing results or proof that an unimplemented candidate layout is integrated into the CLI.

| Requested property | Executed check and observed outcome |
| --- | --- |
| Compare current implementation with alternatives | All five layouts completed all 105 randomized processes. Complete rows, acquisition boundaries and workload checksums matched across layouts. |
| Reflect actual access rather than only a synthetic loop | Actual public tophat calls, including adaptation and sorted cursor reuse, were measured. Public `extractChromatograms` semantic checks passed, and the separate native shared-library tests plus five real extraction CLI/output comparisons passed 14/14. |
| Improve numeric speed where claimed | Across all four large fixture/size cases, factored-SoA dense process-median ranges were entirely below direct-AoS ranges. At 4,224,750 IM points: 4.901-6.995 ms versus 29.759-31.042 ms. |
| Reduce memory where claimed | Across all four large cases, the maximum factored-SoA retained allocation was below the minimum RawTile retained allocation. Large IM: 163.55 versus 550.33 MiB. This excludes retained source and does not claim whole-process reduction. |
| Expose less favorable access costs | Large IM factored AoS beat factored SoA on actual public narrow extraction across all process medians: 115.80-118.11 versus 125.94-174.49 ms. Direct AoS full-row access was faster in every trial. Universal SoA superiority was rejected. |
| Preserve semantics and report uncertainty | Full-field bit-exact reconstruction and per-result extractor comparisons passed. All timing samples and broad noisy trial ranges were retained, with no outlier deletion. |
| Reproduce the reported result | The full 105-process summary was regenerated after completion and matched the committed JSON byte-for-byte. All 65 published table cells were independently matched to that result. |

This closes the **requested comparison experiment** with observed improvements and observed counterexamples. It does not close a production-layout migration: none was requested or implemented, and no full future-engine acceptance path exists yet for these alternative layouts. The lack of that adapter is a scope limit, not a successful full-engine validation or an externally blocked test hidden as a pass.

## Limits and next decision boundary

This measures raw-point layouts only. It does not settle hypothesis/evidence/family/score layouts, byte-budget admission, streaming lifetime, multithreading, large diaPASEF-frame shapes, larger target libraries, full-run throughput or end-to-end scientific results. Million-point inputs repeat whole real acquisitions without inventing new peak distributions. The target set has only 34 queries, so the balance of adaptation versus extraction work can differ for a real library.

There is not yet a production consumer of the v0 records. Consequently the measured real legacy access pattern is a constraint, not a complete trace of the future engine. A future IM scan view could share unit/presence metadata, avoid per-point optional objects and avoid repacking. That could improve both memory and extraction, but is **unmeasured**, not a claim made by this report.

The experiment supports a concrete decision now: separate acquisition metadata from numeric points, preserve full semantics, and provide a simple bulk view. Choose the final IM hot layout when that real adapter exists. Do not add a generic layout framework or rewrite unrelated contracts to chase these microbenchmark ratios.
