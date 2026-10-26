# WP-02: Byte budget and bounded queues

- **Track:** data handling.
- **Status:** private primitives and local acceptance checks implemented and verified 26 October 2026. Keep this package open and stacked pending adapter integration and combined-design validation. **Owner / reviewer:** ownership is tracked in GitHub.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10664](https://github.com/OpenMS/OpenMS/issues/10664).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).
- **Implementation:** `feature/opendia-wp02-resource-budget`, based directly on WP-01 `3cd24a25026615b26671a46596b3eeae7de89986`. Primitives: `6621ef44ab`. Record ownership and integration checks: `338ba2c569`. Backpressure-test observation fix: `4b2c8e1d33`. No merge into `develop` and no physical-layout or scientific-semantics freeze.

## Objective

Bound in-flight work with one byte budget and predictable progress rather than a multi-pool resource framework.

## Scope

- Implement byte reservations, release/cancellation and bounded queues using the private v0 interfaces.
- Define a fixed resource-acquisition order and bounded downstream drain capacity before admitting producers.
- Account live decoder, index, evidence, model and writer buffers as their adapters become available, with a measured runtime/native-memory reserve.
- Expose reserved/committed bytes and queue high-water marks separately from externally measured RSS/cgroup peaks.
- Split, spill or fail oversized work explicitly. Never truncate scientific populations to meet memory limits.
- Defer general disk quotas and multi-volume admission. WP-08 owns v1 free-space checks and safe writes.

## Dependencies and integration policy

- **Start / local readiness:** WP-01 only. Branches and future PRs remain stacked on their prerequisites.
- Resource hardening is not a blocker for WP-03/04/05/08 adapters or their early existing-engine comparisons.
- Attach adapters as they become available. Whole-process budget proof belongs to final integration, not this primitive-level milestone.
- The 26 October owner decision supersedes earlier independent local-merge recommendations. Do not merge work-package PRs into `develop` or close their issues until the combined refactor demonstrates correctness, speed and memory efficiency and the owner approves integration.

## Implemented private surface

- `ByteBudget`: one configured total, fixed native allowance, minimum admission, checked arithmetic, stop-aware waiting and permanent cancellation.
- `ByteReservation`: move-only pending bytes. Admit a complete peak job once, then split its reservation without a second blocking admission.
- `MemoryLease`: move-only committed charge, released automatically, with shared accounting state that can outlive the budget object.
- `BoundedQueue<T>`: preallocated fixed-capacity ring, move-only payload support, graceful close/drain, per-call stop and immediate cancellation. Queue capacity counts items, not payload bytes.
- Optional unique charges attach to `RawTile`, `HypothesisPage` and `EvidenceBlock` shared immutable storage. Copies share one charge. Storage destroys its payload before releasing accounting, including validation failures and the last remaining copy.
- `requiredBytes()` checks retained vector capacity and string capacity, with conservative object/allocation allowances. An undersized nonempty charge is rejected. Validation, decoder/model workspace and by-value accessor copies require separate admission or native headroom.

These remain private, uninstalled v0 contracts. There is no new CLI option, allocator interceptor, resource-pool abstraction, layout rewrite or dependency beyond standard-library threading.

## Fixed acquisition protocol

1. Configure native headroom and reserve fixed decoder/index/model workspace, bounded writer/drain space and queue-node storage **before producers start**. Reject configurations that cannot admit a useful complete job after this startup reservation.
2. Before materialization, atomically reserve the whole next job's peak input, output and scratch bytes.
3. Split the admitted reservation into component charges. Declare charges before their allocations so they are destroyed last. Payload ownership carries its charge through queues and copies.
4. Do not request more blocking admission while holding that job. Consumers reuse pre-admitted drain workspace. Queue flow is acyclic and forward only, so downstream completion releases capacity for upstream admission.
5. Stop producer admission and wake all blocked calls on failure/cancellation. Discard queued payloads outside the queue mutex, join all users, then destroy queues and the budget. Live charges are never revoked early.

This acquisition order is a caller contract, not runtime deadlock detection. Requests outside the configured minimum/usable bounds fail immediately with `InvalidValue`, including after cancellation. Ordinary full admission waits or returns `nullopt` through `tryReserve()`. Stopped/cancelled admission returns `nullopt`. Scheduling fairness and admission deadlines are not promised. The regression tests have explicit deadlines.

## Local readiness checks and observations

All five local checks are satisfied for the private implementation, not for a full rewritten engine:

- [x] Exact byte boundaries, checked overflow including `SIZE_MAX`, invalid/below-minimum configuration and requests, move/split/commit, exceptions and cancellation release.
- [x] A synchronized stalled-writer producer/consumer test terminates within its fixed timeout without exceeding either queue's capacity.
- [x] Contending admissions follow the fixed whole-job protocol; per-call stop and permanent cancellation wake waiters without leaked reservations.
- [x] Oversized admission fails explicitly and leaves the input population unchanged. No spill implementation or scientific truncation is introduced.
- [x] Native tests report reservation/commit/queue peaks, with external process RSS measured separately and no accounting-as-RSS claim.

| Requirement / boundary | Observed result on the final implementation |
|---|---|
| Arithmetic/admission/ownership | Always-active `OpenDIAResourceBudget_test` passes exact/minimum/overflow, partial/whole split, zero/over-commit, replacement, exception, token-outlives-budget, stopped/global cancellation and 8-thread contention checks. |
| Queue lifecycle/backpressure | `OpenDIABoundedQueue_test` passes wraparound/order, move-only ownership, failed-push retention, full/empty stop/close/cancel, reentrant cancellation destructors and contending producers behind a gated consumer. |
| Real data and immutable record integration | `OpenDIAResourceInterop_test` uses the actual `MzMLFile` loader, actual record constructors, public `ChromatogramExtractorAlgorithm` and a bounded writer. It preserves every coordinate/value through the queues and reopens the generated TSV to compare every row. Non-IM: **513/513**, IM: **32,750/32,750**, FAIMS: **44,354/44,354**. The extraction kernel is also checked against a direct numerical sum, without substituting a copied implementation. |
| Stalled drain and cancellation | Both queues reach peak **2** under the declared capacity **2**, with a simultaneously blocked byte admission. Releasing the writer drains the complete population. Cancelling the blocked IM pipeline writes **0** rows and releases all record/job charges. After fixed fixture storage is destroyed, reserved and committed bytes are both **0**. |
| Long/shared payloads and invalid records | Capacity and long native-ID/provenance strings are charged. Copies retain one charge. All three storage types reject a one-byte-undercharged exact moved payload. Malformed/nonfinite/count-capped records release charges while preserving the caller's original population. |
| Actual native project / packaging | Debug GNU **13.3.0**, existing native SDK, `STL_DEBUG=OFF`: all **84 selected native OpenDIA and installed-consumer CTests** pass, including **43 installed-consumer fixture steps**, core-only/CLI layering, no-Boost public headers and relocated consumers. All **6 nested installed-SDK consumers** pass. This is actual project compilation and CLI output regression, not a stub library. |
| Standalone / NDEBUG | Full standalone SDK build passes **8/8** tests. A separate Debug build with explicit `-DNDEBUG` on the real implementation and tests passes **8/8**. Budget/queue tests each pass **50 consecutive runs**. After the follow-up observation correction below, real-resource integration passes **50 consecutive native and 50 consecutive NDEBUG runs**. |
| Sanitizers | Actual private archive and three resource tests instrumented with ASan/UBSan pass **3/3**, with leak detection and no diagnostics/suppressions. Shared OpenMS/core/external dependencies are not transitively instrumented. No ThreadSanitizer result is claimed. |

### Follow-up replay and test-race correction

A fresh native replay on **26 October 2026 at 02:17 UTC** passed **83/84**, exposing an intermittent failure in the integration test rather than accepting an isolated successful retry. The global OpenMS termination message showed the last deliberately caught oversized-admission exception. Reporting the active exception instead identified `Did not exercise full queue and blocked byte admission`, reproduced on the **ninth** diagnostic-only execution.

The test had already witnessed the full gated writer queue and a blocked admission, but then reread the transient waiter count. The extractor could release a raw job between those observations. Commit **`4b2c8e1d33`** latches the witnessed condition and reports active exceptions. It changes only the test, with no timeout extension, removed assertion, retry-to-pass behavior or production accounting/queue change.

The corrected source was rebuilt against the actual native project and SDK. Observed results: **84/84 native/installed-consumer tests**, **8/8 standalone**, **8/8 explicit NDEBUG**, **50 consecutive native and 50 consecutive NDEBUG real-data pipeline executions**, and **3/3 ASan/UBSan resource tests**, all passed. Verbose observations again preserved **513**, **32,750** and **44,354** rows, reached both queue caps of **2**, stayed within **68,961,200** accounted bytes and released all charges. Cancellation wrote **0** rows. These supersede the earlier five-run pipeline stress evidence.

The same follow-up passed **21/21 layout semantic tests** and **14/14 actual native extraction class/CLI/FuzzyDiff tests**. The [layout comparison](../experiments/layout/RESULTS.md) still supports metadata factoring and reports SoA counterexamples, not a universal layout winner.

Actual CLI output and installed-package checks establish compatibility. The bounded integration still uses a fixture-specific TSV writer and declared decoder allowance. No production streaming decoder, full CLI allocation enforcement or full-engine speed/memory improvement is established by these checks.

### External memory is a separate measurement

The earlier native resource-test process on `338ba2c569`, measured with `/usr/bin/time -v`, exited **0** with peak RSS **79,960 KiB**, about **78.1 MiB**. This is a recorded fixture measurement, not a remeasurement after the test-only observation correction. Successful cases each report:

- pending-reservation high-water: **67,178,672 bytes**, including the initial fixed-workspace admission;
- committed high-water: **68,652,976 bytes**;
- simultaneous reserved-plus-committed high-water: **68,961,200 bytes**, equal to the configured usable cap;
- declared native allowance: **134,217,728 bytes** (128 MiB);
- configured total including that allowance: **203,178,928 bytes**;
- raw/writer queue high-water: **2/2 items**.

The reserved and committed high-water values occur at potentially different instants and must not be added. The cancelled case's committed peak is scheduling-dependent, while its accounting cap, queue bounds, zero writes and final zero outstanding charges are asserted.

A separate real `OpenDIA --help` startup process exits **0** and peaks at **60,888 KiB**. The 128 MiB fixture native allowance exceeds these observed startup/test process peaks, but is **not** a derived universal reserve for models, third-party decoders, children or representative searches. The test's fixed 64 MiB decoder/reference allowance is declared before loading and retained until all fixture storage is destroyed. It is not an instrumented bounded mzML decoder.

These are fixture measurements, not speedup, layout adoption, whole-engine RSS enforcement or the final 4 GiB/2 GiB acceptance result.

## Reproduction

Run from the implementation checkout. The existing out-of-tree native build retains its historical WP-01 directory name while compiling the stacked WP-02 branch:

```bash
cmake --build ../OpenMS-wp01-native-build --target OpenDIA OpenDIAEntry_test OpenDIARecords_test OpenDIAQuery_test OpenDIAFamilies_test OpenDIARecordsInterop_test OpenDIAResourceBudget_test OpenDIABoundedQueue_test OpenDIAResourceInterop_test -j3
ctest --test-dir ../OpenMS-wp01-native-build -R 'OpenDIA|TestExternalCode' --output-on-failure -j3
ctest --test-dir ../OpenMS-wp01-native-build/installed-consumer/build --output-on-failure
ctest --test-dir ../OpenMS-wp01-native-build -R 'OpenDIAResource|OpenDIABoundedQueue' -V

cmake -S src/topp/OpenDIA -B ../OpenMS-wp02-build -DCMAKE_BUILD_TYPE=Debug -DOpenMS_DIR="$PWD/../OpenMS-wp01-native-build" -DCMAKE_CXX_COMPILER=g++-13
cmake --build ../OpenMS-wp02-build -j2
ctest --test-dir ../OpenMS-wp02-build --output-on-failure
cmake -S src/topp/OpenDIA -B ../OpenMS-wp02-ndebug-build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS=-DNDEBUG -DOpenMS_DIR="$PWD/../OpenMS-wp01-native-build" -DCMAKE_CXX_COMPILER=g++-13
cmake --build ../OpenMS-wp02-ndebug-build -j2
ctest --test-dir ../OpenMS-wp02-ndebug-build --output-on-failure
ctest --test-dir ../OpenMS-wp02-ndebug-build -R 'ResourceBudget|BoundedQueue' --repeat until-fail:50 --output-on-failure
ctest --test-dir ../OpenMS-wp01-native-build -R ResourceInterop --repeat until-fail:50 --output-on-failure
ctest --test-dir ../OpenMS-wp02-ndebug-build -R ResourceInterop --repeat until-fail:50 --output-on-failure

/usr/bin/time -v ../OpenMS-wp01-native-build/src/tests/class_tests/bin/OpenDIAResourceInterop_test
/usr/bin/time -v ../OpenMS-wp01-native-build/bin/OpenDIA --help
```

The real-resource test removes only its own uniquely named temporary TSV. It does not alter source fixtures or frozen WP-00 artifacts.

## Remaining integration work

- The [WP-03 raw adapter](03-raw-tiles.md#small-fixture-process-measurement) now applies private scan/copy charges in the actual CLI. Nine fresh PASEF fixture runs retained matching output. Median peak RSS was 80,884 KiB for the normal reader, 82,012 KiB at the 65 MiB tiled setting and 81,372 KiB at 256 MiB. A lower private cache cap did not lower process peak RSS on this small fixture.
- A separate local heap profile identifies native metadata, controlled-vocabulary and Xerces decode allocations outside private enforcement. These observations do not validate the fixed 64 MiB native allowance for representative runs. Measure and bound those paths before a whole-process claim.
- The [full TimSim follow-up](03-raw-tiles.md#measured-reader-cost-and-copy-admission) found 14 fitting copies rejected by maximum-size admission. The adapter now reserves the known copy bound before allocation. Actual pinned-copy pressure still fails within the cap. A new 1,024 MiB run completes with matching full outputs, but takes 1,040.29 seconds at 1,830,824 KiB RSS. At 128 MiB, decoder admission still fails when views remain retained. These results do not validate the native allowance, a process RSS limit or an efficiency gain.
- Existing v0 callers may still omit accounting to preserve the compatibility path. The current production OpenDIA CLI is **not** globally governed by this budget.
- WP-03/04/05/08 must supply conservative pre-materialization decoder/index/model/output charges and apply the admission protocol to their actual buffers. A test-only fixed allowance is not that adapter integration.
- WP-09 must connect failure propagation and cancellation through the real raw-to-result workflow. WP-10 must measure representative full-process/child RSS, throughput, disk and scientific equivalence before approving the combined design.
- Windows/macOS, cgroup/child peaks and whole-stack sanitizer coverage have not been exercised. `clang-format` is unavailable in this environment; local whitespace checks pass.

**Current dependency increment:** the first WP-03 indexed mzML adapter is tested on this stacked foundation. Full acquisition and resource checks remain open. WP-04 and WP-08 can proceed alongside it. Keep every package exploratory until the combined-design acceptance decision.
