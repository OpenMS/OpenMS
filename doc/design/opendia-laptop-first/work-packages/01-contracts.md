# WP-01: Contracts and scientific semantics

- **Suggested issue title:** [OpenDIA][WP-01] Define scientific and bounded-output contracts
- **Parent epic:** not created
- **Issue:** not created
- **Repository source:** `doc/design/opendia-laptop-first/work-packages/01-contracts.md`

Repository context: `doc/design/opendia-laptop-first/README.md` and `doc/design/opendia-laptop-first/architecture.md`.

- **Status:** proposed, not implemented.
- **Owner:** unassigned. Suggested expertise: C++ architecture and DIA scientific semantics.
- **Reviewer:** unassigned, with producer/consumer owner sign-off.
- **Dependencies:** WP-00 for fixtures and baseline policy.

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

## Code placement

- Provisional until implementation: generic packed records, IDs and views belong in core `DATASTRUCTURES` (serialization in `FORMAT`), with matching header/source paths.
- OpenDIA candidate semantics and schema/profile policy belong in NEW proposed `src/openms/include/OpenMS/ANALYSIS/OPENDIA/` and `src/openms/source/ANALYSIS/OPENDIA/`, neither currently present. This plan creates no directories or code.
- Core must not depend on OPENDIA. Existing shared `OPENSWATH` components stay shared, and `src/topp/OpenDIA.cpp` remains a thin CLI.

## Output-size contract

- Default results-only schema contains compact identification, localization, quantification and ambiguity provenance, not whole-library traces or transition cubes. Existing full-OSW compatibility is an opt-in export only.
- Own a versioned byte/cardinality ledger: result rows, ambiguity links, distinct dictionary entries, sparse incidence edges and shared coordinate/evidence references, with encoded widths and variable-field lengths. Never duplicate sequence/protein strings or evidence arrays per alternative.
- Require profile manifests to declare row counts, actual encoded bytes, mandatory fields and completeness. Growth follows reported results plus required provenance edges, not library size multiplied by scans/transitions.
- Optional rescore profile declares and preserves all prefilter target/decoy/alternative inputs required by each supported rescoring method. Diagnostics are opt-in, explicitly selected and byte capped.
- Working inference scratch, final deliverables and reusable caches have separate budgets plus a per-volume combined peak. Include WAL/journals, indices, atomic staging and simultaneous old/new compaction copies.
- Preserve the full inference population until global confidence completes. Quotas never authorize silent dropping or evidence pruning. If mandatory output cannot fit, fail explicitly without a complete-looking truncated result.
- Compression or Parquet alone does not satisfy this contract. Any CLI budget names remain proposals, not implemented options.

## PR-sized implementation steps

1. Add a short reviewed contract specification and identity/units tables, with input/output examples from WP-00.
2. Introduce minimal owning records and const views with byte-accounting hooks. Register new public headers and required tests/bindings if exposed.
3. Add policy keys, version validation, stable serialization examples and cancellation/error contracts.
4. Exercise one producer-to-consumer chain with real fixture data. Obtain sign-off from raw, extraction, scoring and storage owners before freezing version 1.

## Acceptance checks

- [ ] Review header/source placement and build dependencies: reusable core has no OPENDIA dependency, shared OPENSWATH stays shared, and the CLI contains only thin orchestration.

- [ ] Round-trip each profile with exact row/edge/dictionary counts and encoded-byte accounting, retaining ambiguity provenance without duplicated evidence.
- [ ] Grow absent hypotheses on a fixed result fixture and show default output has no library-by-scan/transition growth; validate declared rescore inputs separately.
- [ ] Exhaust mandatory-output quota and reject publication as complete; test separate and combined-volume accounting including staging/compaction.

- [ ] Unit tests cover query equality and inequality for every scientific field, overflow/malformed offsets, missing versus zero, and invalid lifetime use with appropriate sanitizers.
- [ ] Round trips preserve exact IDs, coordinate references, masks and unresolved status. Stable IDs do not depend on thread scheduling or tile size.
- [ ] Calibration updates invalidate dependent query/evidence results. Unknown units/schema versions fail explicitly.
- [ ] A real fixture crosses the implemented producer/consumer boundaries without rebuilding full string-keyed maps or materializing the whole library.
- [ ] After WP-09, the public workflow demonstrates the same contract invariants end to end. A synthetic transport test alone does not close that gate.

## Risks and handoff

An over-general object graph can recreate the memory problem. Prefer small packed contracts with cold metadata separated. Provide headers/specification, serialization examples, ownership diagrams and test fixtures to WP-02, WP-03, WP-04 and WP-08.
