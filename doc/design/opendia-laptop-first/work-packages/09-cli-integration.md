# WP-09: Opt-in public CLI and end-to-end integration

- Status: proposed / not implemented.
- Reviewer: unassigned.
- Owner: unassigned. Suggested expertise: TOPP integration, CTest, workflow delivery.
- Context: [master plan](../README.md), [architecture](../architecture.md).

## Goal

Deliver an explicit opt-in public OpenDIA path that connects bounded raw access,
mobility-native retrieval, family scoring and global confidence to usable outputs.
Preserve existing outputs and defaults until scientific-quality and resource gates
pass. A working internal harness is not completion of this work package.

## Dependencies and handoff

- Requires [WP-03 raw tiles](03-raw-tiles.md).
- Requires [WP-04 hypothesis store](04-hypothesis-store.md).
- Requires [WP-05 shared extraction](05-shared-extraction.md).
- Requires [WP-06 candidate retrieval](06-candidate-retrieval.md).
- Requires [WP-07 family scoring](07-family-scoring.md).
- Requires [WP-08 storage/confidence](08-storage-confidence.md).
- Uses acceptance definitions in [WP-00](00-acceptance-baseline.md).
- Hands measured integrated results to [WP-10](10-performance-release.md).

## Scope and non-goals

- Wire acquisition, paged hypotheses, shared queries, inference and global output.
- Expose explicit engine selection, memory/disk budgets, worker count, working
  storage, diagnostics and approximation policy through documented public options.
- Include real mobility and localization in the first integrated vertical slice.
- Preserve raw scoring access and uncertainty recovery across stage boundaries.
- Maintain supported existing outputs and require explicit capability negotiation.
- Do not silently fall back to reduced evidence, tilewise FDR or forced winners.
- Do not switch defaults, claim release readiness or replace WP-10 benchmarking.

## Deliverables and interface ownership

- Own CLI option names, validation, help text, workflow wiring and exit semantics.
- Own public CTest fixtures and reproducible raw-to-export command recipes.
- Keep proposed option names uncommitted until parameter review with maintainers.
- Record source/binary revision, parameters, input hashes, calibration/model versions,
  storage capability, search coverage and completion status in run provenance.
- Validate unsupported capability combinations before writing final artifacts.
- Propagate reservation, spill, cancellation and disk-limit outcomes across stages.
- Expose reserved versus committed memory and observed process-peak diagnostics.
- Coordinate legacy output mappings with WP-08, including ambiguous localization.
- Keep inference/export runnable through supported public entry points, not only
  private test helpers. Publish exact commands when executable options exist.

## Small PR slices

1. Add opt-in parameter/help contracts, preflight checks and negative public tests.
   Keep current execution and default output behavior unchanged.
2. Connect the first raw mobility-to-family-scoring vertical slice, including
   missing-MS1/anchor and uncertain-IM cases plus full raw evidence revisits.
3. Connect working storage, global confidence and existing export adapters.
   Add actual public CLI round trips for localization and ambiguous/coexisting forms.
4. Add overflow/cancellation regressions, resource measurements and reproducible
   acceptance scripts. Hand measured artifacts and open gates to WP-10.

## Requirement-to-acceptance checks

- Public workflow: build affected targets locally from one recorded source revision,
  then run the actual opt-in OpenDIA executable from raw input to final exports.
  CTest must assert nonempty scientific results, not merely successful process exit.
- Native mobility: real acquisition fixtures exercise RT-IM overlap, calibration
  offsets/uncertainty, missing scans, tile boundaries and discrete FAIMS geometry.
  Include non-IM compatibility without reducing the mobility fixture to one axis.
- Family science: known-site standards, unresolved/coexisting mixtures, absent
  hypotheses and no-MS1/no-anchor cases retain expected evidence and ambiguity.
  Measure identification, localization and quantitative error separately.
- Storage independence: repeat supported backend paths through public inference
  and export. Verify localization inputs/results survive and unsupported requests
  fail clearly before final output creation, without a scientific downgrade.
- Global confidence: compare integrated results with WP-08 global-oracle fixtures
  and held-out calibration gates, not tilewise q-values or target/decoy symmetry only.
- Resource/scheduling: run CPU-only 2/4 GiB settings at 2/4 workers, varied tiles,
  dense families and slow writers. Record process/cgroup peak, swaps, temporary
  disk, explicit spill/failure and stable candidate/score/quantity reconciliation.
- Compatibility: run affected existing TOPP output comparisons and opt-in outputs.
  Preserve supported field semantics and default behavior, including failure modes.
- Reproducibility: save exact commands, hashes, logs and comparison tolerances.
  Time cold construction and warm reuse end to end at matched search/quality scope.
  Passing tiny fixtures or extraction timing cannot establish greenfield speed gains.

## Risks and completion handoff

- Older binaries versus newer references are not current-source acceptance.
  Resolve build prerequisites and output regressions before claiming validation.
- Baseline mobility/IPF observations do not prove a joint new-engine workflow.
- Hand WP-10 reproducible commands, quality/resource matrices, compatibility results
  and explicit failed/unmeasured gates. Keep the engine opt-in until release gates pass.
