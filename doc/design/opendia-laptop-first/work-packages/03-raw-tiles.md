# WP-03: Mobility-native bounded raw evidence tiles

Status: proposed/not implemented
Owner: unassigned
Reviewer: unassigned
Suggested expertise: mass-spectrometry acquisition formats, ion mobility, packed storage and I/O.
Index: [Laptop-first design](../README.md)

## Goal

Provide byte-bounded acquisition tiles that preserve the original measurement geometry.
Serve local evidence queries without whole-run retention or a dense RT/mz/mobility cube.
Non-IM, continuous IM, diaPASEF and FAIMS must be explicit supported or rejected modes.

## Dependencies

- Requires [WP-01 contracts](01-contracts.md) and [WP-02 budgets](02-resource-budget.md).
- Feeds [WP-05 extraction](05-shared-extraction.md) and [WP-06 retrieval](06-candidate-retrieval.md).
- Public acceptance is integrated by [WP-09](09-cli-integration.md).
- Use [WP-00 baseline](00-acceptance-baseline.md) fixtures and measurement provenance.

## Scope and non-goals

- Implement bounded decoding, optional indexed disk tiles, range access and shared cache leases.
- Preserve exact times, m/z, intensity, mobility units and calibration/coordinate references.
- Preserve acquisition validity and actual frame/scan isolation membership.
- Do not implement extraction aggregation, calibration fitting, candidate gates or scoring.
- Coarse bins are lookup accelerators only, never replacement measurement coordinates.
- Do not silently turn FAIMS compensation voltage into a continuous mobility axis.

## Deliverables and interface ownership

- Own the raw-store adapter, tile reader/writer and cache implementation under WP-01 schemas.
- Address tiles by run, acquisition lane, scan/frame interval and optional mobility interval.
- Pack numeric arrays with offsets by frame and mobility-scan range, plus exact-time references.
- Keep independent MS1/MS2 axes, including multiple MS1 scans within a nominal cycle.
- Represent absent mobility explicitly for non-IM runs, not as an observed zero coordinate.
- Retain separate FAIMS voltage strata in lane identity and every range query.
- Record diaPASEF isolation geometry per frame/scan, including moving or overlapping windows.
- Expose acquisition-validity masks distinguishing no acquisition from acquired zero signal.
- Preserve separated mobility components for downstream local patches and mobilogram views.
- Validate required acquisition metadata and reject unsupported or ambiguous geometry explicitly.
- Define core/halo ownership, exact boundary coordinates and overlap query semantics.
- A halo does not itself establish equivalent peak picking or complete candidate evidence.
- Provide lifetime-safe immutable views and bounded raw-region access for later scoring.
- Cache entries carry raw content identity, schema and coordinate/calibration provenance.
- Invalidate incompatible entries and detect corrupt or incomplete on-disk tiles on reopen.
- Account decoder scratch, compressed bytes, decoded bytes, pinned views and mmap residency.
- Use WP-02 admission before decoding and bounded spill/index construction, including disk quota.
- Oversized frames split only at valid coordinate boundaries or fail with a resource-limit result.

## Small PR slices

1. Add packed tile/axis/validity representation and round-trip tests for ordinary DIA and IM.
2. Add frame-dependent PASEF geometry and discrete FAIMS lanes with negative metadata tests.
3. Add bounded decoder/cache and indexed disk access with cancellation and corruption checks.
4. Wire real acquisition fixtures into WP-05/WP-09 and compare direct versus cached evidence.

## Acceptance checks

- Unit: round trips preserve exact coordinates, offsets, acquisition masks and lane identities.
- Unit: missing scans differ from acquired empty scans and zero-intensity observations.
- Unit: overlapping PASEF windows, moving boundaries, separate FAIMS voltages and absent IM.
- Unit: interval overlap includes long queries whose centres lie outside the requested tile.
- Unit: cache eviction cannot invalidate a leased view, and cancellation releases pinned bytes.
- Unit: oversized frame, truncated cache, disk exhaustion and missing metadata fail explicitly.
- Integration: rebuild coherently and run the actual public pipeline through WP-09 after integration.
- Exercise real mobility-resolved acquisition with overlapping components and tile-edge signals.
- Include non-IM, PASEF and FAIMS inputs, missing scans and multiple MS1 acquisitions.
- Compare direct decode and disk-cache paths for evidence, candidate sets, scores and quantities.
- Vary tile size and 2/4 workers, checking deterministic ownership and no truncated/duplicate evidence.
- Reuse baseline families `TOPP_OpenSwathWorkflow_17`, `_17_cache` and `_23` as starting fixtures.
- Add positive FAIMS and missing-acquisition fixtures rather than claiming those names cover them.
- Measure whole-process/cgroup peaks, bytes decoded/revisited, cache hits and temporary disk use.
- Proposed 2/4 GiB targets remain unproved until measured across preparation through export.
- Store commands, build IDs and fresh producer/output checks in the WP-00 acceptance ledger.

## Failure risks

- Averaged cycle times or mobility collapse can destroy resolving evidence despite smaller storage.
- A static rectangular window model can assign diaPASEF evidence to the wrong hypothesis.
- A nominally bounded cache can still pin unlimited pages through consumer views or mmap.
- Unsupported vendor metadata must not silently become ordinary DIA or measured zero values.

## Handoff

Give WP-05/WP-06 the tile query API, fixture manifests and coordinate/validity invariants.
Document decoder-specific untracked memory and minimum frame allocations for WP-02.
Keep public-pipeline acceptance open until fresh end-to-end outputs pass, not merely round trips.
