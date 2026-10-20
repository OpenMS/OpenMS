# WP-06: Mobility-native candidate retrieval

- Status: proposed / not implemented.
- Reviewer: unassigned.
- Owner: unassigned. Suggested expertise: DIA search, calibration, indexed retrieval.
- Context: [master plan](../README.md), [architecture](../architecture.md).

## Goal

Retrieve high-recall local candidate families in RT and mobility without requiring
a global feature map, detectable MS1, or a strong shared fragment anchor.
Preserve recovery paths so an early ranking cannot irreversibly choose the form
that later identification and localization are allowed to consider.

## Dependencies and handoff

- Requires [WP-03 raw tiles](03-raw-tiles.md) for acquisition-valid raw queries.
- Requires [WP-04 hypothesis store](04-hypothesis-store.md) for paged families.
- Requires [WP-05 shared extraction](05-shared-extraction.md) for evidence views.
- Inherits contracts and reservations through these dependencies.
- Supplies candidates to [WP-07 family scoring](07-family-scoring.md).
- Public integration belongs to [WP-09](09-cli-integration.md).

## Scope and non-goals

- Combine observed-evidence proposals with targeted hypothesis recovery queries.
- Use calibrated mass, RT and IM bounds with explicit uncertainty provenance.
- Broaden or recover when calibration is weak, missing or out of domain.
- Retain RT-IM patches, separated components and exact acquisition references.
- Respect frame/scan isolation geometry and discrete FAIMS strata.
- Represent absent mobility explicitly for non-IM data using the same contract.
- Do not require an MS1 observation or a universally shared anchor fragment.
- Do not implement final confidence estimation, exporters or model training here.
- Do not treat top-k, beam limits or candidate truncation as lossless budgeting.

## Deliverables and interface ownership

- Extend [WP-01 base contracts](01-contracts.md) with retrieval request/results,
  reviewed with WP-03 through WP-07, rather than independently redefining them.
- Own run mass/RT/IM calibration fitting, bounded anchor selection and uncertainty
  validation. Anchors aid calibration but are not mandatory candidate evidence.
- Freeze/version calibration snapshots and invalidate dependent queries/results
  on updates. Missing or unreliable anchors/calibration require broad fallback.
- Request: tile/acquisition references, calibration version and uncertainty,
  hypothesis scope, query bounds and a reservation-backed work allowance.
- Result: stable family/form IDs, local support, retrieval reasons, alternatives,
  null/decoy membership, recovery state and explicit completeness metadata.
- Keep raw references usable by scoring for full-resolution additional queries.
- Own candidate completeness domains and selection semantics, paged iteration and
  reconciliation of broad peaks across tiles, including support and competing forms.
  A fixed halo or matching IDs alone does not establish complete peak evidence.
- Emit fan-out, rejected/recovered counts and approximation provenance by reason.
- Share repeated evidence through WP-05 query identity rather than precursor ID.
- Split, spill or slow oversized regions without silently deleting alternatives.
- If exhaustive processing cannot fit, return an explicit resource-limit result.

## Small PR slices

1. Define request/result contracts and fixtures for missing MS1, missing anchors,
   uncertain IM and competing sites. Establish an exhaustive small-case oracle.
2. Implement bounded anchor selection, run mass/RT/IM fitting, validated uncertainty
   and frozen/versioned snapshots with explicit invalidation and broad fallback.
3. Add geometry-aware lookup and bounded iteration, including non-IM/FAIMS
   contracts and broad-peak completeness/reconciliation.
4. Add uncertainty expansion and targeted recovery for weak observed proposals.
   Expose rejection reasons and approximation state in machine-readable output.
5. Add shared-evidence scheduling, overflow recovery and recall regression tests.
   Wire an integration fixture for WP-09 without switching production defaults.

## Requirement-to-acceptance checks

- High recall: compare retrieved IDs against exhaustive small-case enumeration.
  Include low abundance, modification classes, absent hypotheses and ambiguity.
  Record pre-FDR recall by stratum with thresholds agreed in WP-00 before tuning.
- Uncertainty fallback: perturb calibration and remove IM predictions. Verify
  broad/recovery paths retain truth rather than reporting a confident empty set.
- Calibration ownership: validate held-out mass/RT/IM residual coverage, bounded
  anchor selection and sparse/no-anchor fallback. Snapshot changes must invalidate
  dependent queries/results while frozen snapshots reproduce the same candidates.
- Native mobility: overlapping RT signals separated in IM remain distinct through
  retrieval, including missing scans, isolation changes and tile boundaries.
- No mandatory anchor/MS1: recover known families with either or both absent.
- Safe selection: apply equivalent target/decoy logic, preserve null hypotheses,
  then use held-out/entrapment checks. Symmetry alone is not proof of valid FDR.
- Scheduling: reconcile candidate sets across tile sizes and 2/4 workers, with
  deterministic IDs and no duplicate boundary hypotheses or unexplained losses.
- Completeness: peaks wider than configured halos retain full support and competing
  forms versus an untiled oracle. Check declared domains, selection and recovery,
  not only ID equality, including peaks spanning more than two adjacent tiles.
- Bounded work: dense regions and large families trigger bounded paging/spill
  or explicit failure. Measure reserved, committed and process peak memory.
- Public acceptance after WP-09: execute the actual opt-in OpenDIA CLI from raw
  mobility input through exported results, inspect retained truth and recovery
  diagnostics, and repeat at 2/4 GiB. Unit fixtures alone cannot close this gate.

## Risks and completion handoff

- Recovery can increase disk/time substantially. Publish coverage and cost, not
  an unsupported claim that narrowing is automatically fast or statistically safe.
- Weak discriminating evidence must remain accessible to WP-07 after retrieval.
- Hand off fixtures, oracle comparisons, calibration assumptions, query contracts
  and unresolved recall strata to WP-07/WP-09 and confidence owners in WP-08.
- Proposed gates are not measured greenfield gains. No default replacement occurs
  before integrated scientific-quality and resource acceptance.
