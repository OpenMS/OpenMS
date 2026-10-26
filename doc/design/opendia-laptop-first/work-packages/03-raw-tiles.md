# WP-03: Bounded raw evidence tiles

- **Track:** data handling.
- **Status:** first indexed mzML adapter tested. Full acquisition and resource gates remain open.
- **Branch:** `feature/opendia-wp03-raw-tiles`, stacked on WP-02 `4b2c8e1d33`.
- **Implementation:** [3a460113a3](https://github.com/OpenMS/OpenMS/commit/3a460113a391a63178d0a6da7da32a24ef63e2c5).
- **IM CLI checks:** [1c46042570](https://github.com/OpenMS/OpenMS/commit/1c46042570cf83a21d92c4968b5b935d5bee5bb4).
- **Query-scale check:** [57381c0213](https://github.com/OpenMS/OpenMS/commit/57381c0213c86e08a6c5ded49f64bf3024ed1d5d).
- **Parent:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10665](https://github.com/OpenMS/OpenMS/issues/10665).
- **Policies:** [architecture](../architecture.md), [master plan](../README.md).

## Objective

Read raw evidence without keeping all run peaks in memory. Preserve acquisition geometry and observable evidence. Do not merge into `develop` before the combined design passes its correctness, speed and memory checks.

## Current design

- Reuse the indexed mzML reader and its SWATH-map checks. Do not add a disk-cache format.
- Use one scan per tile. Store metadata once per tile and numeric data in columns.
- Keep MS1 and MS2 maps separate. Preserve native RT, m/z, IM units, frame/scan metadata and isolation bounds without rounding.
- Keep FAIMS discrete. Do not convert absent mobility to zero or absent scans to zero-intensity scans.
- Share the cache between reader clones. Keep byte charges until the last view is released.
- Give the existing engine separate mutable copies. Charge those copies too.
- Add `readOptions=tiled`. Keep the default reader unchanged.
- Fail if a scan is too large or retained views prevent admission. Do not wait while holding views.

## Checked on 26 October 2026

| Requirement | Real path and observed result |
|---|---|
| Preserve scans, arrays and metadata | `OpenDIARawAccess_test`: exact native metadata and numeric arrays for 114/150/19 scans and 513/32,750/44,354 points from non-IM, IM and FAIMS fixtures. |
| Preserve extraction at different limits | The public chromatogram extractor gives exact RT/intensity arrays at two cache limits and 1/16/4,096 queries. Reference signals are nonzero. |
| Keep raw charges independent of query count | For each fixed traversal and limit, raw peak bytes and scan misses are equal at all three query counts. IM uses 23,080 bytes at the tight limit and 962,424 bytes at the large limit, with 150 misses. Output arrays grow from 300 to 1,228,800 elements. This is a coordinate-grid test, not a complete peptide-library test. |
| Preserve interval edges and missing values | Tests cover long isolation overlap, open/closed RT and isolation edges, absent IM, measured zero IM, discrete FAIMS zero and empty scans. |
| Preserve engine results | Direct and tiled CLI runs match four non-IM rows and six PASEF mzML rows, plus all nonempty `FEATURE*`/`SCORE*` tables, at 65/256 MiB declared budgets. Identity fields are exact. Numeric tolerance is relative `1e-6`, absolute `1e-8`. Only generated run/file/feature IDs are normalized. |
| Exercise IM scores | Both PASEF workflows contain finite, nonzero measured IM and IM intensity scores. The check rejects a copied result with those values removed. Equal nonfinite values are also rejected. |
| Keep owners and charges valid | Tests cover eviction, pinned admission failure, array-only/container-only owners, mutation isolation, cancellation and last-owner release. |
| Reject invalid input | Tests reject invalid indices/queries, oversized scans/copies, a missing file index and a removed indexed payload. The real CLI rejects an indexed fixture with three FAIMS CV lanes without process termination. No partial tile is published. |
| Return parallel reader errors | Six actual indexed-reader cases return admission errors from IM calibration, legacy extraction and wave extraction. The CLI oversize test exits with an error, not a termination signal. |

The first adapter passed 98 native/project/installed-SDK checks, the installed OpenDIA consumer, 10 standalone tests, 10 NDEBUG tests and two ASan/UBSan raw tests. The follow-up added eight real CLI checks. All 106 native/project/installed-SDK checks passed. Instrumentation covers the private implementation, not the complete native SDK. These are Linux checks. `clang-format` is not installed. `git diff --check` passed.

On the non-IM CLI fixture, both tiled budgets reached a private accounted peak of 606,128 bytes. The PASEF fixture used 1,230,760 bytes at the larger budget. The 65 MiB case forced cache eviction and kept private charges within its 1 MiB usable cap. Results still matched. At the run log point, 10,920 non-IM bytes and 23,768 IM bytes remained in live engine-held copies. This is not a leak-free whole-engine or process-memory claim. The tests confirm final-owner release in the private adapter.

These tiny fixtures test software behavior, not error-rate calibration. Percolator reports that their training populations are too small. Do not use their scores as scientific acceptance evidence.

## Small-fixture process measurement

Nine fresh complete PASEF mzML CLI processes ran in randomized order, with three samples per mode. Each used the same prepared library, two threads and the native Debug SDK. All exported rows and feature/score tables match, including nonzero IM evidence. GNU `time` measured process peak RSS separately from private cache bytes.

| Reader | Peak RSS KiB, min / median / max | Wall seconds, min / median / max |
|---|---|---|
| Normal | 80,788 / 80,884 / 81,192 | 1.10 / 1.10 / 1.12 |
| Tiled, 65 MiB declared | 81,876 / 82,012 / 82,116 | 1.07 / 1.09 / 1.11 |
| Tiled, 256 MiB declared | 81,256 / 81,372 / 81,528 | 1.08 / 1.10 / 1.11 |

The smaller cache did **not** reduce process peak RSS on this fixture. These short, non-exclusive-host runs do not establish a speed gain. Do not use them to set a general native reserve or to accept the rewritten-engine resource targets.

A separate local heaptrack run retained matching output. Allocation stacks identify native metadata and controlled-vocabulary loading under `loadOnDisc` and Xerces decoding under `RawAccess`. They are outside private charge enforcement. The profiler reports 15.83 MB peak heap, not RSS. Its RSS includes profiler overhead. Its reported process-exit allocations are not a classified leak result. No whole-engine leak claim is made.

The final test commit passed all 106 native/project/installed-SDK checks, 10 standalone checks and 10 NDEBUG checks. The query-scale assertions ran in each build. Local evidence includes all command vectors, nine samples and the separate heap profile under `wp03-20261026/memory-0326`. No frozen input was changed.

## Limits

- The first adapter supports indexed mzML spectrum maps. Split-window input is not supported. The native loader rejects mixed FAIMS lanes rather than merging them.
- Actual WP-00 input checks reject the TimSim `.d.zip` input and Astral `.raw` input. The first adapter cannot yet run the matched WP-00 resource comparison. This is a missing input path, not a failed speed or memory target.
- Private tile/copy charges are not a process RSS limit. Native metadata, indexes, decoding and control storage use a declared 64 MiB allowance. That allowance is not measured or enforced by an allocator.
- A large native scan can allocate memory before the private size check. Later growth of mutable copies and downstream engine workspaces is not intercepted.
- Cache access uses one lock and linear LRU lookup. No whole-engine speed or memory improvement is claimed.
- Oversized scans fail. Split scans only at valid acquisition boundaries if measured input requires this.

## Remaining gates

- [ ] Test native per-frame/scan diaPASEF geometry and lane-aware FAIMS through the complete engine. PASEF mzML parity and mixed-CV rejection do not close this gate.
- [ ] Test complete peptide libraries and representative raw files. The 4,096-query public-extractor check supports fixed-traversal raw-charge independence, but does not bound output arrays or full-engine hypothesis memory.
- [ ] Add supported representative input paths, or pin equivalent indexed mzML conversions for both readers. Keep conversion cost separate. Then measure native decoder/metadata and downstream workspace memory, runtime and peak RSS against WP-00. Do not change the frozen baseline files.
- [ ] Supply the private scan view to WP-05/06. Complete WP-09/10 integration and acceptance checks.

Keep this issue open. Keep WP-01, WP-02 and WP-03 stacked. Store planning documents on the planning branch only.
