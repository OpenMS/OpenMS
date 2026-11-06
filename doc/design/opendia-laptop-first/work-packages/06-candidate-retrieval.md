# WP-06: Mobility-native candidate retrieval

- **Track:** science.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10668](https://github.com/OpenMS/OpenMS/issues/10668).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Retrieve complete local candidate families with high recall, preserving recovery when observed evidence or calibration is weak.

## Grouped peak proposal, 5 November 2026

The current workflow's per-transition peak picks and later consensus/overlap reconciliation motivate a group-level candidate proposal before expensive downstream scoring. A robust group trace may nominate shared RT boundaries, but it must not become a hard filter based on a fixed transition majority: weak or absent fragments, transition-specific delays, interference and partially coeluting groups can otherwise lose true candidates. Treat group support as a cheap proposal/priority signal, retain a broad/recovery route for weak or discordant evidence, and measure candidate recall against an exhaustive small-case oracle before using it to reject work. Keep per-transition evidence available for scoring and explanation. The first implementation is an OpenDIA-scoped opt-in path; other OpenSWATH workflows retain legacy behavior.

## Initial grouped-picker benchmark, 6 November 2026

A single matched prepared-PQP run on Astral RAW 01 with 32 threads completed end to end. Relative to the active-RT repeat-2 baseline, full wall time was 13:49.15 versus 13:17.69 (+31.46 s, about 4.0% slower); the candidate RT scoring stage was 9:45 versus 9:22 (+23 s). Peak RSS was essentially unchanged (8,355 versus 8,361 MiB). The export contained 1,511 rows versus 1,501. Matching by precursor ID, 1,487 rows were shared, 14 baseline-only and 24 grouped-picker-only. Several shared-row quantities and RT boundaries changed, so this is not output parity and does not establish improved candidate recall. This is one run per mode, not a robust timing study. Keep `robust_consensus` opt-in, not the OpenDIA default, and do not use it to reject downstream candidates until exhaustive-oracle recall and stratified scientific validation pass. This bridge did not demonstrate a speed improvement on this workflow.

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
- [ ] Group-proposed candidates preserve recall for missing/weak fragments, one-transition interference, delayed transitions, overlapping true peaks and non-identical RT sampling; fallback/recovery reasons are recorded.
- [ ] Removing MS1/shared anchors or perturbing calibration retains truth through recorded broad/recovery queries.
- [ ] Peaks spanning more than two tiles preserve full support and competing forms versus an untiled oracle across tile sizes/workers.
- [ ] IM-separated components, moving isolation windows and discrete FAIMS fixtures retain expected acquisition membership and candidate identities.
- [ ] Target/decoy/null counts and selection reasons are reproducible, and oversized-region tests page/spill or fail without silent truncation.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
