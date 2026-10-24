# WP-04: Bounded hypothesis preparation and indexed family store

- **Track:** library.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10666](https://github.com/OpenMS/OpenMS/issues/10666).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Prepare and query complete hypothesis populations without materializing the whole library or localization graph in memory.

## Scope

- Stream imported libraries and FASTA digestion, modifications, decoys, prediction batches and index construction.
- Keep stable family/form/precursor/fragment IDs, target/decoy provenance and cold protein metadata with paged lookup.
- Store shared fragments and sparse incidence once, without pairwise UIS expansion or hidden top-k/modification truncation.
- Preserve recovery access for missing or uncertain RT/IM predictions and broad interval-overlap queries.
- Plug a paged-store adapter opt-in into the existing engine before depending on new family science.
- Record cold preparation versus reuse costs, model startup, page peaks, store cardinality and actual disk use.

## Dependencies

- **Start / local merge:** WP-01 only.
- WP-02 hardening is not a blocker. WP-05 tile/store adapters integrate as available, not as prerequisites here.
- Feed WP-06 complete paged families, and provide stable references for WP-08.

## Local merge criteria

- [ ] Streamed digestion/modification/decoy enumeration exactly matches an exhaustive small fixture, including protein mappings and alternative forms.
- [ ] Changing page/batch/worker sizes preserves all population IDs and counts, including an oversized family processed in multiple pages.
- [ ] Existing-engine adapter on/off runs yield identical candidate identities and tolerance-equal scores/quantities for imported and library-free fixtures.
- [ ] Missing predictions and interval queries whose centres lie outside the region retain the expected recovery hypotheses.
- [ ] Interrupted or incompatible stores are rejected on reopen, and a scale test records bounded live pages plus preparation time and bytes without hidden truncation.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
