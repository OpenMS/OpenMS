# WP-03: Mobility-native bounded raw evidence tiles

- **Track:** data handling.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10665](https://github.com/OpenMS/OpenMS/issues/10665).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Replace whole-run raw retention with bounded access while preserving observable evidence and acquisition geometry.

## Scope

- Add a tool-private raw adapter, packed tiles, range queries and leased cache views.
- Preserve exact RT/mz/IM coordinates, independent MS1/MS2 axes and per-frame/scan diaPASEF isolation geometry.
- Keep FAIMS discrete, absent mobility explicit, and missing acquisition distinct from measured zero.
- Add an opt-in adapter to the existing engine first, comparing direct and tiled paths before new retrieval/scoring.
- Cache size follows raw observations and indexes, not hypothesis count. Oversized frames split at valid boundaries or fail explicitly.

## Dependencies

- **Start / local merge:** WP-01 only.
- WP-02 hardening is not a prerequisite. Integrate its admission hooks as available.
- Supply tile access to WP-05/06 without waiting for their completed science or final WP-09 integration.

## Local merge criteria

- [ ] Round-trip fixtures preserve exact coordinates, acquisition masks and lane identity for non-IM, diaPASEF and discrete FAIMS.
- [ ] Long interval-overlap queries and tile-edge signals match direct access without lost or duplicated observations.
- [ ] Existing-engine adapter on/off runs produce identical identities and tolerance-equal evidence, scores and quantities on pinned fixtures.
- [ ] Eviction and cancellation tests preserve leased views and release all pinned buffers. Corrupt/truncated caches fail explicitly.
- [ ] Varying tile sizes and hypothesis counts leaves unique raw-observation counts unchanged, with recorded cache bytes and bounded live tiles.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
