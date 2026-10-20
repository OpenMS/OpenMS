# WP-06: Mobility-native candidate retrieval

- **Track:** science.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10668](https://github.com/OpenMS/OpenMS/issues/10668).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Retrieve complete local candidate families with high recall, preserving recovery when observed evidence or calibration is weak.

## Scope

- Combine observed proposals with targeted recovery over paged families and full-resolution raw evidence.
- Consume WP-05's existing mass/RT/IM calibration adapter and immutable snapshots, rather than owning a prerequisite calibration rewrite.
- Define uncertainty expansion, retrieval reasons, completeness domains and cross-tile reconciliation.
- Preserve competing/null hypotheses without requiring detectable MS1 or a shared anchor.
- Keep RT/IM-separated components and acquisition eligibility intact. Broaden uncertain queries instead of confidently excluding them.
- Restrict new family retrieval science to single-run v1. Existing multi-input behavior is not redesigned here.

## Dependencies

- **Start / local merge:** WP-03, WP-04 and WP-05.
- Supply candidates to final WP-07 integration. WP-07's early localization experiment does not depend on this package.

## Local merge criteria

- [ ] Candidate IDs and pre-FDR recall match an exhaustive small-case oracle within predeclared thresholds, stratified by abundance and ambiguity.
- [ ] Removing MS1/shared anchors or perturbing calibration retains truth through recorded broad/recovery queries.
- [ ] Peaks spanning more than two tiles preserve full support and competing forms versus an untiled oracle across tile sizes/workers.
- [ ] IM-separated components, moving isolation windows and discrete FAIMS fixtures retain expected acquisition membership and candidate identities.
- [ ] Target/decoy/null counts and selection reasons are reproducible, and oversized-region tests page/spill or fail without silent truncation.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
