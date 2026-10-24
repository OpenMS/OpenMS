# WP-02: Byte budget and bounded queues

- **Track:** data handling.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10664](https://github.com/OpenMS/OpenMS/issues/10664).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Bound in-flight work with one byte budget and predictable progress rather than a multi-pool resource framework.

## Scope

- Implement byte reservations, release/cancellation and bounded queues using the private v0 interfaces.
- Define a fixed resource-acquisition order and bounded downstream drain capacity before admitting producers.
- Account live decoder, evidence, model and writer buffers, with a measured runtime/native-memory reserve.
- Expose reserved/committed bytes and queue high-water marks separately from externally measured RSS/cgroup peaks.
- Split, spill or fail oversized work explicitly. Never truncate scientific populations to meet memory limits.
- Defer general disk quotas and multi-volume admission. WP-08 owns v1 free-space checks and safe writes.

## Dependencies

- **Start / local merge:** WP-01 only.
- Resource hardening is not a blocker for WP-03/04/05/08 adapters or their early existing-engine comparisons.
- Attach adapters as they become available. Whole-process budget proof belongs to final integration, not each component's merge.

## Local merge criteria

- [ ] Unit tests cover exact byte-budget boundaries, overflow, below-minimum requests and release after exceptions/cancellation.
- [ ] A producer/consumer stress test with a stalled writer terminates within a fixed timeout and never exceeds configured queue bounds.
- [ ] Contending tasks follow the documented fixed acquisition order, and cancellation wakes blocked waiters without leaked reservations.
- [ ] Oversized work tests exercise split/spill or explicit failure and assert unchanged input population counts.
- [ ] A local harness reports reservation and queue peaks, plus external process memory, without presenting accounting totals as an RSS guarantee.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
