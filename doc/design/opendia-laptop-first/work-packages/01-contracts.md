# WP-01: Contracts and scientific semantics

[Plan index](../README.md) | [Architecture](../architecture.md)

- **Status:** proposed, not implemented.
- **Owner:** unassigned. Suggested expertise: C++ architecture and DIA scientific semantics.
- **Reviewer:** unassigned, with producer/consumer owner sign-off.
- **Dependencies:** [WP-00](00-acceptance-baseline.md) for fixtures and baseline policy.

## Goal

Define the minimal stable boundaries that let raw access, hypothesis preparation, extraction, candidate generation, scoring and storage be implemented independently without sacrificing memory ownership or scientific meaning.

## Scope and non-goals

Own shared records, IDs, views, ownership, policy versions and error/completeness states. Do not dictate an optimized kernel, final storage backend or new statistical model. Names below are provisional until reviewed.

## Deliverables and owned interfaces

| Contract | Required semantics |
|---|---|
| `RawTile` and raw-region access | Run/frame/scan identity, exact RT/mz/IM coordinates and units, isolation geometry, MS level, FAIMS stratum, validity and leased ownership |
| `HypothesisStore` records | Stable family/form/precursor/fragment IDs, target/decoy provenance, ragged incidence, cold metadata and paged access |
| `ExtractionQuery` | Complete equality/hash key, bounds and edge policies, aggregation, calibration snapshot, isotope and acquisition selection |
| `EvidenceBlock` | Packed ragged/sparse support and masks, distinct missing/zero states, shared immutable axes and lifetime-safe views |
| `CandidateFamily` | Region, alternatives/nulls, shared/discriminating evidence, boundary ownership and completeness domain |
| `ScoreBatch` | Stable result identities, separate identification/localization/quantification state, raw references and provenance |

Agree numeric precision and overflow limits for offsets/IDs. Document what happens when a view outlives a tile, when cancellation occurs, when memory reservation fails and when input is malformed. Serializers must reject incompatible schema or policy versions rather than reinterpret them.

Define calibration as an immutable snapshot within a pass. Cache keys include its identity. Explicitly specify fallback or failure when calibration is uncertain or unavailable.

## Scientific decisions to record

- Tolerance units and full/half widths, open/closed boundaries, sum/max policy, duplicate centroids and coordinate precision.
- Acquisition membership for overlapping/PASEF windows and discrete FAIMS, with absent mobility allowed.
- Exact query sharing versus merely overlapping evidence. Shared observations are not statistically independent duplicates.
- Candidate completeness and cross-tile reconciliation. A halo and unique apex owner alone are insufficient.
- Reference compatibility versus intentional algorithm changes. Legacy whole-assay picker totals/overlap suppression must be preserved only when claiming reference semantics.
- How unresolved localizations, coexisting forms, absent MS1 and incomplete searches appear in records and outputs.

## PR-sized implementation steps

1. Add a short reviewed contract specification and identity/units tables, with input/output examples from WP-00.
2. Introduce minimal owning records and const views with byte-accounting hooks. Register new public headers and required tests/bindings if exposed.
3. Add policy keys, version validation, stable serialization examples and cancellation/error contracts.
4. Exercise one producer-to-consumer chain with real fixture data. Obtain sign-off from raw, extraction, scoring and storage owners before freezing version 1.

## Acceptance checks

- Unit tests cover query equality and inequality for every scientific field, overflow/malformed offsets, missing versus zero, and invalid lifetime use with appropriate sanitizers.
- Round trips preserve exact IDs, coordinate references, masks and unresolved status. Stable IDs do not depend on thread scheduling or tile size.
- Calibration updates invalidate dependent query/evidence results. Unknown units/schema versions fail explicitly.
- A real fixture crosses the implemented producer/consumer boundaries without rebuilding full string-keyed maps or materializing the whole library.
- After [WP-09](09-cli-integration.md), the public workflow demonstrates the same contract invariants end to end. A synthetic transport test alone does not close that gate.

## Risks and handoff

An over-general object graph can recreate the memory problem. Prefer small packed contracts with cold metadata separated. Provide headers/specification, serialization examples, ownership diagrams and test fixtures to [WP-02](02-resource-budget.md), [WP-03](03-raw-tiles.md), [WP-04](04-hypothesis-store.md) and [WP-08](08-storage-confidence.md).
