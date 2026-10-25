# WP-01: Contracts and scientific semantics

- **Track:** data handling.
- **Status:** in progress, first private-engine increment started 25 October 2026. **Owner / reviewer:** ownership is tracked in GitHub.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10663](https://github.com/OpenMS/OpenMS/issues/10663).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Provide small, explicitly unstable v0 interfaces that unblock independent adapters before scientific contracts are frozen.

## Scope

- Define raw tiles, paged hypotheses, complete extraction keys, sparse evidence, candidate families and score records.
- Include IM coordinates, absent mobility, family alternatives and nulls from day one, without requiring finished IM/family science for the first runnable engine.
- Keep OpenDIA code tool-private under `src/topp/OpenDIA/{OpenDIA.cpp,include,source,test,CMakeLists.txt}`.
- Use one private static/OBJECT library shared by the tool and tests, without installed OpenDIA headers or a new core API.
- Specify ownership, units, stable identities, calibration snapshot identity and schema-version rejection, not a general capability-manifest framework.

## Implementation roadmap

Deliver small, reviewable commits in the order below. The first increment is in
progress on `feature/opendia-wp01-private-engine`, based on clean upstream
`3e47bc7e7468656d0d53237f11fc4f049d278efe`. Keep implementation commits separate
from this scratch-only planning branch. Do not include the optional WP-00 profiler.

### 1. Private engine boundary (in progress)

- Relocate the current tool into `src/topp/OpenDIA/` without changing its CLI,
  parameter defaults, scientific algorithms or output policy.
- Put the existing orchestration behind one private STATIC/OBJECT target. Both
  the thin executable and tool-private tests must use that same implementation.
- Retain standalone tool builds against the OpenMS CMake package. Do not install
  or export an OpenDIA library/header or add a core/pyOpenMS API.
- Check the shared entry point, parameter/INI compatibility and invalid CLI paths.
  Build out of tree and rerun the existing public OpenDIA TOPP producers,
  SQLite/Parquet comparisons, predicted mode and failure-path tests.

### 2. Coordinates, identities and ownership (pending)

- Add minimal v0 raw tile, hypothesis-page and sparse evidence value records.
  Document units and exact coordinate references, stable run/precursor/family IDs,
  calibration identity and explicit schema version.
- Represent absent mobility, measured zero versus missing evidence, and discrete
  FAIMS strata without substituting lookup bins for measurement coordinates.
- Start with bounded owning containers and explicit checked views/leases. Avoid
  raw pointers that can outlive their owner and avoid a general storage framework.
- Tests transfer records between an actual producer and consumer, retain owned
  data across producer destruction, and reject invalid offsets, dimensions,
  schema versions and invalid views. Round trips here need not invent a wire format.

### 3. Complete extraction query identity (pending)

- Define a value key covering run/lane and acquisition selection, MS level and
  isolation geometry, RT/mz/IM bounds and units, discrete FAIMS selection,
  tolerance, edge/aggregation policy, isotope policy and calibration snapshot.
- Use exact equality on canonical scientific values. A hash is an index hint,
  not proof that two queries can share evidence. Reject nonfinite/invalid inputs.
- Vary every field independently in tests. Force hash collisions and prove that
  unequal queries do not share. Calibration changes must invalidate equivalence.
- Do not implement a cache or extraction algorithm in this increment.

### 4. Families, scores and v0 handoff (pending)

- Add minimal candidate-family and score-batch records with stable references,
  unresolved/coexisting alternatives and null hypotheses. Keep identification,
  localization and quantity distinct. State candidate completeness explicitly.
- Retained scores/evidence must not depend on a destroyed raw tile or local worker.
  Do not freeze new family partitioning, null/decoy science or confidence semantics.
- Exercise the complete v0 producer/consumer path with non-IM, IM and discrete
  FAIMS examples, absent versus zero support, alternatives and malformed records.
- Publish a short contract document with field units, ownership, validators and
  schema rejection. Mark scientific freeze pending the WP-07 experiment review.
  Hand v0 to WP-02/03/04/05/08 and early WP-09 without waiting for that freeze.

## Requirement-to-check map

| Requirement | Planned acceptance check |
|---|---|
| One real private implementation, no public API | Native target graph and link commands show executable/test linking the same target. Build-tree and standalone tool builds work. Installation/export files contain no OpenDIA interface. |
| Behavior-preserving first move | Existing TOPP workflow tests pass with the relocated binary. Parameter/INI output matches the baseline apart from declared executable/build provenance fields. |
| Coordinate and missingness semantics | Producer/consumer record round trips preserve IDs, exact RT/mz/IM coordinates, absent/zero masks and discrete FAIMS strata. |
| No false query sharing | Per-field key mutations and a constant-hasher container retain unequal scientific queries separately. |
| Checked ownership/schema/bounds | Invalid lifetime/view access, offsets, nonfinite coordinates and unknown versions fail explicitly. Owned output survives producer/tile destruction. |
| Family ambiguity preserved | Round trips retain all declared alternatives/nulls and separate identification/localization/quantification fields. Freeze remains pending WP-07, not a local closure blocker. |

## Dependencies

- **Start / local merge:** none, including no WP-00 dependency.
- Initial v0 contracts can close before science review. Contract freeze is a later milestone, only after the independent WP-07 early localization experiment is reviewed.
- WP-02/03/04/05/08 and early WP-09 consume v0 without waiting for that freeze.

## Local merge criteria

- [ ] An out-of-tree build links the tool and private tests against the same static/OBJECT implementation, with no installed OpenDIA interface.
- [ ] Round trips preserve IDs, exact coordinate references, missing versus zero, discrete FAIMS strata and unresolved family alternatives.
- [ ] Query-key tests vary every scientific field, including calibration identity, and reject false sharing under forced hash collisions.
- [ ] Tests reject malformed offsets, unknown schema versions and invalid ownership/lifetime use.
- [ ] A producer/consumer fixture passes using v0 records, and the contract document marks freeze pending WP-07 experiment review rather than blocking local closure.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
