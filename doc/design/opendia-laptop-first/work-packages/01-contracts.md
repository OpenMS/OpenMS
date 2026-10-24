# WP-01: Contracts and scientific semantics

- **Track:** data handling.
- **Status:** proposed. **Owner / reviewer:** unassigned.
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
