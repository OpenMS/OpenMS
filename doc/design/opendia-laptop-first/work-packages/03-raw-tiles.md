# WP-03: Bounded raw evidence tiles

- **Track:** data handling.
- **Status:** first indexed mzML adapter tested. Full acquisition and resource gates remain open.
- **Branch:** `feature/opendia-wp03-raw-tiles`, stacked on WP-02 `4b2c8e1d33`.
- **Implementation:** [3a460113a3](https://github.com/OpenMS/OpenMS/commit/3a460113a391a63178d0a6da7da32a24ef63e2c5).
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
| Preserve extraction at different limits | The public chromatogram extractor gives exact RT/intensity arrays at two cache limits and 1/16 queries. Reference signals are nonzero. |
| Preserve interval edges and missing values | Tests cover long isolation overlap, open/closed RT and isolation edges, absent IM, measured zero IM, discrete FAIMS zero and empty scans. |
| Preserve engine results | Direct and tiled CLI runs match four exported rows and all nonempty `FEATURE*`/`SCORE*` tables at 65/256 MiB declared budgets. Identity fields are exact. Numeric tolerance is relative `1e-6`, absolute `1e-8`. Only generated run/file/feature IDs are normalized. |
| Keep owners and charges valid | Tests cover eviction, pinned admission failure, array-only/container-only owners, mutation isolation, cancellation and last-owner release. |
| Reject invalid input | Tests reject invalid indices/queries, oversized scans/copies, a missing file index and a removed indexed payload. No partial tile is published. |
| Return parallel reader errors | Six actual indexed-reader cases return admission errors from IM calibration, legacy extraction and wave extraction. The CLI oversize test exits with an error, not a termination signal. |

Local builds passed 98 native/project/installed-SDK checks, the installed OpenDIA consumer, 10 standalone tests and 10 NDEBUG tests. Two ASan/UBSan raw tests passed. Instrumentation covers the private implementation, not the complete native SDK. These are Linux checks. `clang-format` is not installed. `git diff --check` passed.

On the CLI fixture, both tiled budgets reached a private accounted peak of 606,128 bytes. At the run log point, 10,920 bytes remained in live engine-held copies. This is not a leak-free whole-engine or process-memory claim. The tests confirm final-owner release in the private adapter.

## Limits

- The first adapter supports indexed mzML spectrum maps. Split-window input is not supported. The native loader rejects mixed FAIMS lanes rather than merging them.
- Private tile/copy charges are not a process RSS limit. Native metadata, indexes, decoding and control storage use a declared 64 MiB allowance. That allowance is not measured or enforced by an allocator.
- A large native scan can allocate memory before the private size check. Later growth of mutable copies and downstream engine workspaces is not intercepted.
- Cache access uses one lock and linear LRU lookup. No whole-engine speed or memory improvement is claimed.
- Oversized scans fail. Split scans only at valid acquisition boundaries if measured input requires this.

## Remaining gates

- [ ] Test native per-frame/scan diaPASEF geometry and lane-aware FAIMS through the complete engine. The current IM/FAIMS raw tests do not close this gate.
- [ ] Test library-scale hypothesis counts and representative raw files. Confirm that raw storage follows observations, not hypothesis count.
- [ ] Measure and bound native decoder/metadata and downstream workspace memory. Compare full-workflow runtime and peak RSS with WP-00.
- [ ] Supply the private scan view to WP-05/06. Complete WP-09/10 integration and acceptance checks.

Keep this issue open. Keep WP-01, WP-02 and WP-03 stacked. Store planning documents on the planning branch only.
