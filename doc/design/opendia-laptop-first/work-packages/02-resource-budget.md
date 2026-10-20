# WP-02: Global resource budgets and deadlock-free admission

Status: proposed/not implemented
Owner: unassigned
Reviewer: unassigned
Suggested expertise: C++ concurrency, allocators, Linux memory accounting and external-memory pipelines.
Index: [Laptop-first design](../README.md)

## Goal

Make resource admission an enforceable whole-pipeline contract, not a per-worker cache hint.
Bound preparation, decoding, extraction, scoring and output without dropping scientific work.
The proposed CPU-only laptop acceptance settings are 2/4 GiB and 2/4 workers.
These are unproved engineering targets, not measured capabilities or throughput promises.

## Dependencies

- Requires [WP-01 contracts](01-contracts.md), including resource-limit and completeness states.
- Baseline measurements come transitively from [WP-00](00-acceptance-baseline.md).
- Supplies [WP-03](03-raw-tiles.md), [WP-04](04-hypothesis-store.md),
  [WP-05](05-shared-extraction.md) and [WP-08](08-storage-confidence.md).
- Public end-to-end acceptance closes with [WP-09](09-cli-integration.md).

## Scope and non-goals

- One shared controller covers raw/index/evidence/candidate/model/writer allocations.
- Include decoder scratch, compressed and decoded copies, queued results and encoder buffers.
- Track preparation and spill disk quotas alongside managed memory reservations.
- Do not implement raw readers, prediction algorithms, scientific pruning or global inference here.
- Do not equate managed bytes with RSS or promise that mmap enforces a resident cap.

## Deliverables and interface ownership

- Own the implementation of WP-01 resource leases, admission outcomes and usage snapshots.
- Expose acquire/release/cancel semantics, requested bytes, reserved capacity and committed bytes.
- Report pool minimums, allocation granularity, in-flight peaks and observed resident usage separately.
- Enforce `raw + index + evidence + candidates + models + writer <= budget - runtime_reserve`.
- Include allocator overhead and native/model uncertainty in a documented safety reserve.
- Constrain admission by explicit configuration and applicable cgroup limits, not host free RAM alone.
- Expose process RSS and cgroup peak/limit measurements, including mapped/file-backed residency.
- Grow committed blocks on demand rather than preallocating the maximum on tiny inputs.
- Reject a request below feasible minimum capacity before allocating, with the minimum explained.
- Never round an allocation above its pool limit or use an unaccounted heap fallback.
- Define one atomic multi-pool reservation operation for each executable work unit.
- A waiting task must not hold partial leases needed by another waiting task to make progress.
- Reserve bounded downstream drain capacity before admitting producers, including ordered output.
- If estimates grow, release/replan, split or spill without changing statistical context.
- An indivisible oversized job returns an explicit resource-limit/incomplete result, never success.
- Cancellation releases leases and wakes waiters, with bounded cleanup and no lost ownership.
- Publish byte/time/disk telemetry so consumers can attribute stalls and cold preparation cost.

## Small PR slices

1. Implement lease accounting, granularity/minimum checks and deterministic unit tests.
2. Add atomic admission, cancellation and producer/drain progress tests with bounded queues.
3. Add spill quota accounting and process/cgroup telemetry with unsupported-platform reporting.
4. Connect stage adapters and public stress fixtures through WP-09 without scientific policy changes.

## Acceptance checks

- Unit: exact budget edges, overflow-safe arithmetic, double release and exception cleanup.
- Unit: tiny inputs commit only needed blocks, and below-minimum requests fail before allocation.
- Unit: contention, nested requests, cancellation and stalled writers terminate without deadlock.
- Unit: oversized estimates trigger split/spill/failure and never silently increase capacity.
- Integration: use a coherent local build and the actual public OpenDIA pipeline after WP-09 wiring.
- Run preparation through export at proposed 2/4 GiB and 2/4 workers with nonempty outputs.
- Measure process/cgroup peak usage externally as well as internal accounting, with no benchmark swap.
- Grow absent hypotheses, localization alternatives and raw density independently on a fixed run.
- Stall the real writer and exhaust temporary disk, checking bounded memory and explicit failure.
- Compare candidate identities, scores and quantities across budgets/workers under declared tolerances.
- Record source/build IDs, commands, stage peaks, disk peaks and outcomes in the WP-00 ledger.
- A unit-only harness or compliant managed-byte counter cannot close whole-process acceptance.
- If native memory exceeds the target, reduce admission/reserve or fail the gate, not relabel RSS.

## Failure risks

- Hidden model or decoder allocations can defeat the reserve and must be measured in integration.
- Ordered completion can retain later results indefinitely unless drain capacity is protected.
- Huge families may require disk/time beyond the requested budget even with bounded resident state.
- Cgroup page-cache charging differs from process RSS, so report both and the measurement scope.

## Handoff

Give downstream owners lease examples, allocation minima, progress invariants and overflow outcomes.
Provide reproducible adversarial tests and telemetry field definitions to WP-09 and WP-10.
Mark implementation acceptance pending until real preparation-to-export resource tests pass.
