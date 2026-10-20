# WP-05: Exact-query shared extraction and calibration adapter

- **Track:** data handling.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10667](https://github.com/OpenMS/OpenMS/issues/10667).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Share identical extraction queries in the existing engine without changing evidence or waiting for new retrieval science.

## Scope

- Own the adapter to existing mass/RT/IM calibration, with immutable per-pass snapshots and cache invalidation on revision.
- Implement complete query equality, collision-safe sharing and immutable sparse evidence with exact axes and validity masks.
- Retain raw-region access, acquisition geometry and complete picker input/state where existing-engine equivalence requires it.
- Plug extraction opt-in into the existing engine and compare sharing on/off before new candidate/family scoring.
- Handle shared tophat first/last-centroid bugs as a standalone core fix on `develop`, with both IM and non-IM tests, not a compatibility mode.
- Keep that mathematical correction separately reviewed and use the same corrected reference for adapter comparisons.

## Dependencies

- **Start / local merge:** WP-01 only. No wait on WP-06 for calibration.
- WP-03 tile and WP-04 store adapters integrate as available. Existing inputs support initial extraction work.
- WP-02 hardening is not a blocker. Feed calibrated evidence and raw revisits to WP-06/07.

## Local merge criteria

- [ ] Hand-computed tests cover first/last centroids, empty/singleton scans, duplicate mz and exact/near boundaries for both IM and non-IM tophat extraction.
- [ ] Every query-key field and forced hash collisions are tested, and changed calibration snapshots invalidate affected evidence.
- [ ] Existing-engine sharing on/off runs retain identical IDs and tolerance-equal evidence, scores and quantities against the corrected reference.
- [ ] Repeated consumers share one payload, while lifetime/cancellation tests release buffers only after the last lease.
- [ ] Tile-edge, non-monotone interval, missing-acquisition and FAIMS fixtures match scalar reference extraction without coordinate collapse.
- [ ] The separate core-fix revision and its IM/non-IM test results are recorded independently from the private adapter change.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
