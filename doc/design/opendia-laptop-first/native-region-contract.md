# Native column-region contract: first implementation slice

This is an OpenDIA-owned point-storage/query contract, not a completed engine, frame decoder, calibration adapter or claim of faster identification. The default CLI still uses the previous access/extraction/scoring path.

## Inspiration and deliberate differences

Inspected `okohlbacher/dia-trace-extract-prototype` at `f1441c3cd83c94249f8fd94a2130c1d7b164be93`, particularly `RunGrid.h`, `Extractor.cpp`, `docs/EXTRACTION.md` and `docs/10-release-report.md`. Borrow the ideas of acquisition-local contiguous columns, shared axes and batched searches, not source code or a new dependency. That prototype reports extraction-only results, not independent full-engine identification/error-control acceptance.

Do not adopt its run-wide pool, fixed-point m/z, cycle-mean time as the only axis, nearest-window-only assignment or mobility collapse as OpenDIA correctness assumptions. The first slice retains exact native double m/z, float intensity and float point mobility, individual acquisition RTs and duplicate points. Optional lossy representations or changed scientific reductions need separate quality validation.

## Concrete boundary

Implementation belongs in `src/topp/OpenDIA/{include,source,test}/OpenDIARawRegion*`.

| Contract | First-slice behavior |
| --- | --- |
| Identity and geometry | Caller supplies run/lane/calibration context and explicit native source indices. Each requested index produces one row in caller order, including repeats, with native ID, exact RT and MS level. Source indices are not invented frame/scan IDs. Full isolation/window geometry and lane selection remain the caller's responsibility. |
| Layout | CSR-like acquisition offsets into contiguous m/z, intensity and point-IM arrays. No per-point strings or spectrum objects, no inheritance from `ISpectrumAccess`, no `SpectrumPtr` payload. Native decoding currently bridges through `IMSExperimentRandomAccess`/`MSSpectrum`; that bridge is replaceable. |
| Precision and missingness | Preserve double m/z and native float values. Point-mobility unit is explicit. Missing point IM is distinct from zero; an unsupported IM-filtered request fails before output mutation. Optional FAIMS zero is distinct from absent FAIMS. Scalar mobility annotations in milliseconds or inverse reduced units are retained separately, never substituted for point IM; unsupported scalar units reject. |
| Admission | Caller bounds acquisition count, total points and native-ID bytes per acquisition. Checked arithmetic admits all retained array capacities before allocation; oversized sources fail without a partial region or silent truncation. Arrays never grow. Fixed allocation is 16 bytes per point, including reserved IM capacity even when absent, plus row/ID/storage costs. Decoder, sorting, allocator and control workspace are not a process-memory guarantee. |
| Ownership | Immutable copied region/scan/range handles retain the same charge. A returned span borrows its owning view's lifetime. Cancellation prevents publication/new admission, not reading already published views. Final owner releases the retained charge. |
| Query | Per-acquisition binary m/z range with explicit open/closed endpoints. Batch sum accepts unsorted physical requests and exact-sized caller output. Validate the whole request batch before writing. No dense query matrix, hidden cache or hot-path allocation. Summation is double in stable source m/z order, not a new max/centroid/interpolation rule. |

Sorted native input must avoid an unnecessary sort/permutation. Unsorted input uses a stable permutation preserving duplicate order and aligned IM. The contract retains only the declared point columns and acquisition identity, not all auxiliary native arrays or full acquisition metadata.

## Current backend ceiling and next integration gate

The existing native Bruker random reader decodes and sorts a full frame for each per-window request under a shared native-handle mutex. Compact output columns alone do not remove those repeated frame decodes or that lock. Do not call the bridge a frame-native backend or infer lower disk I/O from retained bytes.

A full native CPU profile first identified a smaller immediate bottleneck: 75.61% of self user-CPU samples in `SpectrumAccessQuadMZTransforming::getSpectrumById`, dominated by repeated shared-owner array getters inside the calibration loop. `c64e647997` removes that churn without changing calibration arithmetic. This is a generic shared hot-path correction, not evidence that the new region already drives the CLI or that a neural classifier is needed.

Next connect one calibrated extraction/scoring consumer to the region contract, retaining complete wider-IM evidence. Profile that real path before choosing frame/block reuse, fused traversal or a new scorer. The replaceable native producer must eventually reuse a decoded physical region across relevant isolation/IM requests where measurements justify it. Account for index construction and workspace. Do not hide a converted run or truncate candidates. Compatibility changes require the unchanged twelve-output CLI comparison; deliberately new model science follows its independent quality gates.

## Validation scope

Use the actual out-of-tree vendor-enabled project build. Check real repository mzML/IM fixtures and direct native Bruker and sampled Thermo inputs through the producer and query consumer, comparing values/identity and ordinary interior-window sums with the existing native reader and actual extractor. Literal edge tests cover duplicate and endpoint behavior separately from the known legacy tophat edge defects. Check malformed requests, missing/unit-mismatched IM, empty scans, limits, overflow, cancellation, unchanged outputs on failure and retained-view/final-release behavior.

Cold source-open/index construction is not isolated by the current checks. Optional native output includes retained capacity and repeated batch-query diagnostics, not a comparative throughput result. These diagnostics are not an end-to-end engine speedup. Full Thermo science, calibration/scoring integration, automatic low-cap admission, whole-process resource bounds and competitive engine acceptance remain open.

## Observed requirement-to-check mapping, 27 October 2026

Implementation `6495579dde` adds only the three private contract files and their two CMake registrations. The actual vendor-enabled Release build passed six contract/resource tests, independently replayed by the coordinator. After the separate calibration fix, six affected source/extractor/scorer/region tests pass in 3.94 s.

| Requirement / public output | Executed check and observation |
| --- | --- |
| Exact native identity, RT, m/z, intensity and mobility | Real indexed mzML/IM fixtures read through `MzMLFile`/`OnDiscMSExperiment`; every decoded value, identity and row order compared. Native Bruker `.d.zip`: 25 maps / 75 sampled acquisitions. Native Thermo `.raw`: 301 maps / 903 sampled acquisitions. Both region/query paths pass against the actual native source and `ChromatogramExtractorAlgorithm`, without conversion. |
| Bounded retained storage and lifetime | Fixed-capacity tests cover arithmetic overflow, admission refusal, oversized decoded points/IDs and cancellation. Copies of region/scan/range retain immutable values beyond region and budget lifetime; last-owner retained charge returns to zero. Decoder/sort/control workspace is outside that charge. |
| Range and batch semantics | Literal tests cover every endpoint policy, duplicates, stable unsorted source order, repeated source indices, unsorted requests and empty scans. Ordinary interior-window sums also match the real extractor with and without point IM. Missing point IM, mismatched units and an invalid later request reject before any caller output changes. Scalar IM and FAIMS missing/zero remain distinct. |
| Native execution and failure reporting | Both readers enabled in the real build. Initial Thermo call failed runtime discovery and is retained; the retry uses the already configured .NET runtime via `DOTNET_ROOT` and exits zero. The complete sampled Thermo check takes 29.98 s and peaks at 4,353,832 KiB RSS, not a low-memory release claim. |
| Compatibility of the separate shared calibration fix | Quadratic and ppm arithmetic, clone behavior, empty spectra, unchanged source arrays and intensities pass in the actual transformation test. Existing raw mutable-isolation/fallback, extraction, scoring and native failure tests also pass. Full CLI observation is recorded in [WP-05](work-packages/05-shared-extraction.md#native-region-contract-and-calibration-hotspot). |

Local evidence lives under `wp05-data-contract-20261027`: `native-region-acceptance.log`, `native-region-thermo.log`, `calibration-build.log`, `calibration-tests.log` and the frozen/optimized native profile records. The first profile had concurrent-library provenance failure. The second completed science but its DWARF report observer timed out; `recovered-observation.json` separately verifies unchanged frozen library hashes, all twelve outputs, private cap/final release and the recovered flat report, leaving the original failed observer unchanged.
