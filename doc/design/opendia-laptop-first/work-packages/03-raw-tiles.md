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

## Representative input and indexed IM repair

A vendor-enabled OpenMS build converted the complete TimSim `control_01.d.zip` to a separate indexed mzML file. The frozen archive hash did not change. Conversion took 16.33 seconds and used 299,904 KiB peak process RSS. Conversion cost is separate from search cost.

The new file has 734 spectra, 58,884,702 points and 24 isolation windows. A read-only audit matched native frame point sums, RT and isolation bounds. All spectra have standard inverse reduced IM arrays. The matched prepared library has 2,000 targets, 2,000 decoys and 24,000 transitions.

This real input exposed a crash in default transition filtering. The indexed decoder dropped standard IM array names. The full mzML reader did not. [The shared decoder repair](https://github.com/OpenMS/OpenMS/commit/6b4c23d35fa841fac2eb307a4655094df016f129) uses the same CV rule as the full reader and validates drift arrays before filtering. [The parallel filter repair](https://github.com/OpenMS/OpenMS/commit/8e998694e176b522366aa4a76a80b87f6c31229e) returns reader errors to the caller. Scheduling is unchanged.

| Changed path | Observed check |
|---|---|
| Standard IM names, units and values | Real DOM decoder tests pass for standard and non-standard arrays. |
| Missing or unequal drift arrays | Validation tests throw an error instead of reading invalid storage. Empty arrays and closed interval edges pass. |
| Indexed drift filtering | Direct and tiled reads agree on a standard CV array written by the native mzML writer. |
| Parallel reader failures | Twelve actual indexed-reader admission failures return through filtering, calibration and extraction. No retained cache charge remains. |
| CLI and installed package | All 113 affected native/SDK checks pass, plus the installed OpenDIA consumer. |
| Private builds | All 10 standalone Debug and 10 explicit NDEBUG checks pass. |

The earlier stopped normal trial and crashed tiled trial remain failed evidence. Three new attempts used the same full input and library, eight threads and the native Debug SDK. Default transition filtering and auto-IRT stayed enabled. Each service had a 900-second runtime limit and a separate 12 GiB safety limit.

| Reader | Full-input result |
|---|---|
| Normal | GNU `time` reports exit 0, 856.82 seconds and 1,848,496 KiB peak RSS. Output checks pass: 4,688 result rows and 985/985/251 precursor/peptide/protein matrix rows. IM evidence is nonzero. |
| Tiled, 256 MiB declared | The independent 900-second limit stopped the process during extraction. No final output or complete process-RSS record exists. |
| Tiled, 128 MiB declared | Exit 8 during calibration: `Mutable scan admission failed`. This is an error return, not a crash. No final output exists. |

The tool killed the first measurement supervisor at 600 seconds. The normal engine service continued under its own limit. A separate supervisor recovered its complete GNU `time` record and checked its final outputs. The original wrapper exit is unknown. That interruption is retained in the evidence, not reported as a clean benchmark run.

Before/after source, binary, input and library hashes match. The converted input SHA256 is `0e99f830a2b774571ffe3ac8caab4d1e2330a8e3889442ea0954a142df681325`. The library SHA256 is `e62a8484a7b385549f9605fb3a34379f941ff72153218e8a93c9c69daab61173`. Local commands, logs and recovery records are under `wp03-20261026/representative-1312`.

These two tiled attempts did not pass the full-input check. Neither produced complete output. They do not establish a speed gain or RSS reduction.

### Measured reader cost and copy admission

A 300-second CPU sample of the full input collected 19,567 samples with no losses. The indexed DOM decoder had a 49.03% inclusive share. Sampled extraction threads were often asleep on futexes. This does not identify each wait as a cache-lock wait. The shared cache lock still covers decoding and copies. The stopped sample is not a performance result.

GDB captured 14 copy-admission errors. Checked copy bounds were 12,619,207 to 15,999,247 bytes. Free capacity was 19,196,041 to 33,172,568 bytes. Each copy fit, but the code requested the 33,554,432-byte maximum scan size.

[The fix](https://github.com/OpenMS/OpenMS/commit/3a46947b8517207f1d4b4657d54ae823de290396) reserves the known copy bound before allocation. Decoder admission still reserves the maximum scan size. No cap, lock or science setting changed. [The test](https://github.com/OpenMS/OpenMS/commit/84619d94a5dad6423b0b4444dbd6c70864dd85eb) fails on the old code and passes on the fix.

| Requirement | Observed result |
|---|---|
| Admit a fitting copy | The real indexed-reader test returns a detached copy without eviction. |
| Preserve data and storage ownership | Mutation isolation, array-only and container-only ownership, pinned pressure and final charge release pass. |
| Reject actual capacity pressure | Pinned copies reach admission failure without cap overrun or a retained reservation. |
| Preserve supported interfaces | All 113 native/SDK checks, 10 standalone checks and 10 NDEBUG checks pass. |
| Check private memory use | Both raw ASan/UBSan tests pass. Leak detection is off. The core SDK is not instrumented. |

A fresh full-input run used the tested commits, the same full input and library, eight threads and unchanged science settings. Each engine ran in a separate bounded service. No build or profile ran during measurement. Source, binary, SDK, input and library hashes match before and after.

| Reader | Observed result |
|---|---|
| Tiled, 1,024 MiB declared | Engine and wrapper exit 0. Wall time 1,040.29 seconds. Peak process RSS 1,830,824 KiB. All 4,688 result rows, nonzero IM evidence, nonempty feature/score tables and 985/985/251 matrix rows match the completed normal reference. |
| Tiled, 128 MiB declared | Exit 8 after 112.20 seconds during calibration. The fitting-copy fault is no longer the reported error. The decoder reports `Scan admission failed: pinned views or cancellation`. No final output exists. Its partial RSS is not a completed-run comparison. |

The completed tiled run reports 10,771 cache hits, 2,474 misses and 2,143 evictions. Private peak charges are 1,006,624,847 bytes, below the 1,006,632,960-byte usable cap. At the log point, 125,104,072 bytes remain in engine-held copies. Final-owner release is tested separately.

**Full output parity now passes for this one input at 1,024 MiB. Efficiency acceptance does not.** The earlier normal record is 856.82 seconds and 1,848,496 KiB RSS. It is not a clean matched timing trial because its supervisor was interrupted. These single Debug runs show no clear resource benefit. The 128 MiB path still fails. Do not widen its cap or shorten the scientific population to call it a pass.

Local commands and results are under `representative-1312/copy-acceptance-1445`.

### Release completed calibration scans

The next full-input diagnostic found unused MS2 owners during MS1 calibration reads. At the first failed decoder request, 30,264,025 bytes were free. The decoder needed 33,554,432 bytes. The cache was empty and cancellation was off. The diagnostic wrapper exited 0, but the engine exited 8. This is failure evidence, not a completed search.

[The owner-release fix](https://github.com/OpenMS/OpenMS/commit/d7f5972f79fbf423778371036cc28c06edabd81f) drops completed fragment scans before the precursor read. It covers both IM and m/z calibration. Fit inputs, traversal, byte caps and science settings are unchanged.

[The native regression](https://github.com/OpenMS/OpenMS/commit/008696b0f337fb6d9cdac4ddee847f53008dcd4f) runs both real calibration APIs through the indexed reader. Both paths fail on the old core. Both pass after the fix. Fits and positive extraction windows match the unrestricted native reference within `1e-10`. Peak charges equal the 2,512-byte usable cap. Final charges are zero.

All 114 affected native/SDK checks, six nested consumers, 10 standalone checks and 10 explicit NDEBUG checks pass. Both private raw ASan/UBSan tests pass. Leak detection is off, and the core SDK is not instrumented. These checks prove the local lifetime repair, not complete resource acceptance.

A new full-input check used the same input, full library, eight threads and science settings. Source, binary, SDK, input and library hashes match before and after. No build or profile ran during the check.

| Reader | Observed result after owner release |
|---|---|
| Tiled, 1,024 MiB declared | Engine and wrapper exit 0. Wall time 924.74 seconds. Peak process RSS 1,833,572 KiB. All 4,688 result rows, nonzero IM evidence, feature/score tables and 985/985/251 matrix rows match the previous complete output. |
| Tiled, 128 MiB declared | Engine and wrapper exit 8 after 101.62 seconds. Decoder admission still fails during calibration. No final output exists. The partial RSS is not a completed-run comparison. |

The completed run reports 10,760 cache hits, 2,485 misses and 2,159 evictions. Private peak charges are 1,006,631,890 bytes, below the 1,006,632,960-byte usable cap. Engine-held copies retain 125,104,072 bytes at the log point. The local tests check final-owner release separately.

This change removes a proved local lifetime fault and preserves complete outputs. It does not establish lower process RSS or a speed gain. The timings are single Debug observations, not randomized matched trials. Concurrent whole-job admission, native decoder bounds and shared extraction remain open. Do not raise the low-budget cap or reduce the library to call it a pass. Commands and results are under `representative-1312/calibration-acceptance-1530`.

Debug call stacks also include full-array sortedness precondition checks during extraction. Do not use their cost to choose a Release cache policy. Keep the shared lock until optimized, matched measurements support a different policy.

## Optimized full-input comparison

A separate native Release build used source `008696b0f3`. Core and private sources use `-O3 -DNDEBUG`. The executable loads the new Release libraries, not the Debug SDK. All 120 affected native/CLI/SDK checks and six installed consumers pass.

Nine fresh processes used the same complete TimSim input, full library, eight threads and science settings. Each reader occupied each trial position once. Each process had an eight-core CPU quota, a 12 GiB safety limit and a 1,200-second time limit. No build or profile ran during measurement. These are warm-cache trials on a shared host.

| Order | Reader | Wall seconds | Peak RSS KiB | Result |
|---|---|---:|---:|---|
| 1 | Normal | 161.36 | 1,833,760 | Complete, outputs match |
| 2 | Tiled, 256 MiB | 391.16 | 971,932 | Exit 8, no final output. Partial measurement only. |
| 3 | Tiled, 1,024 MiB | 169.75 | 1,728,084 | Complete, outputs match |
| 4 | Tiled, 256 MiB | 386.36 | 956,156 | Complete, outputs match |
| 5 | Tiled, 1,024 MiB | 191.14 | 2,017,308 | Complete, outputs match |
| 6 | Normal | 187.59 | 1,833,752 | Complete, outputs match |
| 7 | Tiled, 1,024 MiB | 165.25 | 1,978,572 | Complete, outputs match |
| 8 | Normal | 208.80 | 1,833,640 | Complete, outputs match |
| 9 | Tiled, 256 MiB | 426.54 | 954,560 | Complete, outputs match |

All eight completed runs preserve 4,688 result rows, all feature/score tables, finite nonzero IM and 985/985/251 matrix rows. Each matches the first optimized normal run. Private peaks stay within 201,326,592 usable bytes at 256 MiB and 1,006,632,960 bytes at 1,024 MiB. Source, executable, loaded dependencies, build flags, commands, input, library, comparator and frozen archive hashes are unchanged.

Trial 2 returns `Scan admission failed: pinned views or cancellation` during extraction. It is an engine failure, not a timeout. The two later successes do not replace it. The predeclared check rejects this matrix and withholds a performance summary. No speed gain or consistent RSS reduction is established at 1,024 MiB. The two completed 256 MiB runs use less RSS but take longer, and that setting does not complete reliably.

**Next:** apply whole-job admission before concurrent decoder/copy work. Then repeat the failed command and the complete matrix without changing the cap or scientific population. Keep the current lock until further measurements justify a change. Native allocation bounds and shared extraction remain separate open work. Local commands and all nine records are under `representative-1312/optimized-setup` and `optimized-acceptance`.

## Limits

- The first adapter supports indexed mzML spectrum maps. Split-window input is not supported. The native loader rejects mixed FAIMS lanes rather than merging them.
- The current WP-03 build has the Thermo and Bruker readers disabled. A separate vendor-enabled OpenMS build can read these formats and write indexed mzML. The tiled adapter still requires indexed mzML. Do not interpret that adapter restriction as a pyOpenMS format restriction.
- Private tile/copy charges are not a process RSS limit. Native metadata, indexes, decoding and control storage use a declared 64 MiB allowance. That allowance is not measured or enforced by an allocator.
- A large native scan can allocate memory before the private size check. Later growth of mutable copies and downstream engine workspaces is not intercepted.
- Cache access uses one lock and linear LRU lookup. No whole-engine speed or memory improvement is claimed.
- Oversized scans fail. Split scans only at valid acquisition boundaries if measured input requires this.

## Remaining gates

- [ ] Test native per-frame/scan diaPASEF geometry and lane-aware FAIMS through the complete engine. PASEF mzML parity and mixed-CV rejection do not close this gate.
- [x] Match complete-library outputs for the converted TimSim control run at the 1,024 MiB setting. This closes one correctness check, not representative resource acceptance.
- [ ] Reduce repeated decoding and retained-copy costs. The checked-copy and calibration-owner repairs preserve complete output parity. The 128 MiB decoder admission still fails. Integrate whole-job admission and shared extraction without changing the scientific population. Query-grid tests do not bound full-engine output or hypothesis memory.
- [ ] Extend pinned input coverage beyond the converted TimSim control run. Test Astral and the remaining cohort. Keep conversion cost separate. Measure native decoder/metadata, downstream workspace memory, runtime and peak RSS. Do not change the frozen WP-00 files.
- [ ] Supply the private scan view to WP-05/06. Complete WP-09/10 integration and acceptance checks.

Keep this issue open. Keep WP-01, WP-02 and WP-03 stacked. Store planning documents on the planning branch only.
