# WP-07: Shared-evidence family scoring and localization

- Status: proposed / not implemented.
- Suggested issue title: `[OpenDIA][WP-07] Shared-evidence family scoring and localization`.
- Parent epic: OpenDIA laptop-first redesign (issue link pending).
- Issue: not created. This document is a posting draft.
- Repository source: `doc/design/opendia-laptop-first/work-packages/07-family-scoring.md`.
- Reviewer: unassigned.
- Owner: unassigned. Suggested expertise: peptidoform inference, DIA/IM scoring.
Repository context: `doc/design/opendia-laptop-first/README.md` and `doc/design/opendia-laptop-first/architecture.md`.

## Goal

Score local families using shared and lazily requested discriminating evidence,
preserving ambiguous or coexisting forms instead of forcing a single localized
winner. Identification, localization and quantification must retain their own
meaning and uncertainty while using native RT-IM evidence from the first slice.

## Dependencies and handoff

- Requires WP-05 shared extraction.
- Requires WP-06 candidate retrieval.
- Consumes their bounded evidence views, complete query identities and raw access.
- Supplies inference records to WP-08 storage/confidence.
- Shares integrated public checks with WP-09.

## Scope and non-goals

- Represent backbone/composition families with explicit site alternatives and
  initial charge/adduct scope. Define any later charge aggregation separately.
- Query shared fragments once and lazily request site-discriminating fragments,
  including informative loss channels and weak distinguishing evidence.
- Use sparse fragment-to-hypothesis incidence and streamed sufficient statistics.
- Preserve RT-IM shape, precursor/fragment agreement and interference information.
- Permit identification without detectable MS1 or a shared anchor.
- Preserve unresolved groups and multiple forms when supported by evidence.
- Quantify explicitly selected RT/IM support without inventing resolved amounts
  for inseparable components. Record ambiguous allocation rather than hiding it.
- Do not implement global q-value reduction or storage-specific inference forks.
- Do not precompute all pairwise signatures or dense forms-by-fragments patches.

## Deliverables and interface ownership

- Extend WP-01 base contracts with scoring inputs/results and
  family-level evidence attribution, rather than independently redefining them.
- Inputs include stable IDs, calibration/acquisition references, sparse incidence,
  optional precursor evidence and a full raw-query capability via WP-05/WP-03.
- Results include identity/localization statistics, competing/coexisting forms,
  ambiguity groups, quantification regions and evidence/provenance references.
- Keep shared observations from being counted as independent replicated evidence.
- Expose missing evidence separately from measured absence and incompatible data.
- Define additional-query requests so scoring can revisit the full raw signal,
  not merely the fragments or compressed summaries selected during retrieval.
- Own bounded per-family accumulation and approximation diagnostics.
- Coordinate serialized inference inputs with WP-08 before fixing schema versions.

## Code placement

- Proposed NEW directories, currently absent: `src/openms/include/OpenMS/ANALYSIS/OPENDIA/`
  and `src/openms/source/ANALYSIS/OPENDIA/` own OpenDIA family workflow,
  attribution/selection policy, ambiguity handling and evidence-query orchestration.
- Reuse generic math in `MATH/STATISTICS` and reusable scoring algorithms in
  `ANALYSIS/OPENSWATH` or appropriate core modules. Generic serialization belongs
  in `FORMAT`/`DATASTRUCTURES`/`SYSTEM`, never dependent on OPENDIA types/policies.

## Output-size contract

- Compact results are the default: stable IDs, dictionaries and sparse shared
  attribution preserve identification, localization and quantitative ambiguity,
  confidence fields and run provenance without a full library copy or dense traces.
- Optional rescore evidence must cover every declared candidate population required
  by supported global rescoring, including rejected targets/decoys and alternative
  forms. Accepted winners or sufficient statistics for one model are not automatically
  sufficient for another. WP-08 declares exactly which future operations are supported.
- Explicit diagnostic traces are selected/capped. Profiles change retention only,
  not evidence requests, hypotheses, scoring or current-run scientific inputs.
- Default storage is lossless. Share evidence rather than silently quantizing,
  pruning or truncating to fit quotas. Streaming/Parquet alone is not a size solution.
- Keep inference scratch until global FDR separate from final artifacts and reusable
  raw/library caches, with explicit quotas and peak per-volume totals covering
  WAL/journal, indexes, staging and compaction coexistence. Mandatory output overflow
  fails before validated final commit, never by discarding scientific alternatives.
- Scratch removal follows global inference, validated commit and declared retention
  needs. Capability manifests distinguish results-only from supported rescore modes
  and require preserved matching raw data/models for referenced future re-extraction.

## Small PR slices

1. Define family scoring contracts and truth fixtures for known sites, ambiguous
   mixtures, absent MS1/anchors and two overlapping mobility components.
2. Implement shared evidence attribution and lazy discriminating queries with
   sparse incidence. Test against explicit exhaustive small-family scoring.
3. Add joint RT-IM/interference scoring and unresolved/coexisting-form outputs,
   including region-based quantification and explicit missing-evidence handling.
4. Add streamed large-family accumulation, spill/revisit paths and integration
   records for WP-08. Document and test any optional beam/top-k approximation.

## Requirement-to-acceptance checks

- [ ] Review includes and build dependencies: family/attribution policy stays in
  OPENDIA, reusable scoring stays in core, and core has no OPENDIA dependency.
- [ ] Round-trip ambiguous/coexisting forms and shared attribution in compact output,
  then rescore the declared full population, including rejected targets/decoys/alts.
- [ ] Compare all retention profiles for identical current-run scores/quantities;
  stress family fan-out and quotas, asserting explicit failure rather than lost forms.
- [ ] Shared evidence: duplicating a hypothesis reference must not increase the
  statistical weight of its shared observation. Compare with a small exact oracle.
- [ ] Lazy evidence: a weak site-discriminating ion changes the supported localization
  when appropriate, even if another form initially ranks higher on shared ions.
- [ ] Full raw access: request an ion/region absent from retrieval summaries and verify
  scoring receives original-resolution signal with correct acquisition eligibility.
- [ ] Native mobility: known overlapping RT/IM components remain distinguishable,
  while unresolved support produces ambiguity rather than a fabricated winner.
- [ ] Optional evidence: test no MS1, no common anchor, full precursor evidence and
  absent IM. Disabling available scientific evidence is not a memory optimization.
- [ ] Calibration: on held-out known-site mixtures, report localization error apart
  from sequence-identification error and quantitative error with WP-08.
  Require predeclared WP-00 thresholds, not agreement with the highest score alone.
- [ ] Approximation: compare any pruning with exhaustive truth, disclose missed forms
  and localization-error changes, and never label truncation as complete inference.
- [ ] Resources: stress large families and slow consumers at 2/4 GiB and 2/4 workers.
  Verify bounded state, reconciled boundary results and explicit overflow outcomes.
- [ ] Public acceptance after WP-09: run actual opt-in OpenDIA raw-to-export CLI on
  mobility-resolved known-site and unresolved mixtures, then rescore stored inputs.
  Compare identities, ambiguity groups, quantities and separate confidence fields.

## Risks and completion handoff

- Exact localization can be expensive even with shared storage. State the search
  coverage and time/disk costs of every approximation or fallback.
- Dependence among fragments can invalidate naive confidence calculations.
- Hand off attribution rules, raw revisit tests, truth fixtures and inference
  inputs to WP-08, plus public integration assertions to WP-09.
- No greenfield accuracy or speed gain is established by these proposed checks.
