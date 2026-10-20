# WP-09: Opt-in public CLI and end-to-end integration

- Status: proposed / not implemented.
- Suggested issue title: `[OpenDIA][WP-09] Opt-in CLI with compact output profiles`.
- Parent epic: OpenDIA laptop-first redesign (issue link pending).
- Issue: not created. This document is a posting draft.
- Repository source: `doc/design/opendia-laptop-first/work-packages/09-cli-integration.md`.
- Reviewer: unassigned.
- Owner: unassigned. Suggested expertise: TOPP integration, CTest, workflow delivery.
Repository context: `doc/design/opendia-laptop-first/README.md` and `doc/design/opendia-laptop-first/architecture.md`.

## Goal

Deliver an explicit opt-in public OpenDIA path that connects bounded raw access,
mobility-native retrieval, family scoring and global confidence to usable outputs.
Preserve existing outputs and defaults until scientific-quality and resource gates
pass. A working internal harness is not completion of this work package.

## Dependencies and handoff

- Requires WP-03 raw tiles.
- Requires WP-04 hypothesis store.
- Requires WP-05 shared extraction.
- Requires WP-06 candidate retrieval.
- Requires WP-07 family scoring.
- Requires WP-08 storage/confidence.
- Uses acceptance definitions in WP-00.
- Hands measured integrated results to WP-10.

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

## Code placement

- Keep `src/topp/OpenDIA.cpp` a thin CLI calling OpenDIA orchestration in proposed
  NEW, currently absent `src/openms/include/OpenMS/ANALYSIS/OPENDIA/` and
  `src/openms/source/ANALYSIS/OPENDIA/`. Add no engine implementation in `topp`.
- Generic math remains in `MATH/STATISTICS`, reusable scoring in `ANALYSIS/OPENSWATH`
  or appropriate core, and generic serialization/sorting in `FORMAT`/`DATASTRUCTURES`/
  `SYSTEM`. No reusable core may depend on OPENDIA types or policies.

## Output-size contract

- For the new opt-in path, expose compact results as the default retained profile,
  preserving identification, localization/quantification ambiguity and run provenance.
  Keep legacy defaults unchanged. Use IDs/dictionaries/sparse shared evidence, with
  no default full-library duplicate or dense traces and no hidden lossy encoding.
- Document optional rescore evidence covering all declared populations required by
  supported global rescoring, including rejected targets/decoys/alternatives, and an
  explicit selected/capped diagnostic traces profile. Profiles never alter run science.
- Expose separate budgets/locations for inference scratch through global FDR, final
  artifacts and reusable raw/library caches, plus peak per-volume totals including
  WAL/journal, indexes, staging and compaction coexistence. Streaming/Parquet is not
  itself a size guarantee. Publish WP-08 footprint audit and WP-00 size-goal checks.
- Preflight optional full OSW compatibility export, including coexistence space.
  Mandatory-output overflow must fail clearly before committed final output, never
  silently dropping hypotheses, truncating or changing scientific inputs to fit.
- Expose capability manifests and actionable errors: results-only cannot promise
  rescore, and raw re-extraction requires preserved matching data/models behind
  references. Document supported public reopen/rescore/export commands and prerequisites.
- Cleanup follows global inference, validated commit and declared retention needs,
  not CLI exit alone. Make failed versus committed artifacts unambiguous to users.

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

- [ ] Review `src/topp/OpenDIA.cpp` for thin option/dispatch code only: engine work
  calls OPENDIA orchestration, with no new engine code in TOPP or core-to-OPENDIA dependency.
- [ ] Through public CLI commands, reopen every profile and compare current-run science,
  ambiguity/provenance and declared rescore populations; reject results-only rescore
  and missing raw/model prerequisites with actionable capability errors.
- [ ] Exercise separate quotas, shared-volume peaks, capped diagnostics and optional
  OSW preflight; assert no committed partial artifact or silent scientific reduction.
- [ ] Public workflow: build affected targets locally from one recorded source revision,
  then run the actual opt-in OpenDIA executable from raw input to final exports.
  CTest must assert nonempty scientific results, not merely successful process exit.
- [ ] Native mobility: real acquisition fixtures exercise RT-IM overlap, calibration
  offsets/uncertainty, missing scans, tile boundaries and discrete FAIMS geometry.
  Include non-IM compatibility without reducing the mobility fixture to one axis.
- [ ] Family science: known-site standards, unresolved/coexisting mixtures, absent
  hypotheses and no-MS1/no-anchor cases retain expected evidence and ambiguity.
  Measure identification, localization and quantitative error separately.
- [ ] Storage independence: repeat supported backend paths through public inference
  and export. Verify localization inputs/results survive and unsupported requests
  fail clearly before final output creation, without a scientific downgrade.
- [ ] Global confidence: compare integrated results with WP-08 global-oracle fixtures
  and held-out calibration gates, not tilewise q-values or target/decoy symmetry only.
- [ ] Resource/scheduling: run CPU-only 2/4 GiB settings at 2/4 workers, varied tiles,
  dense families and slow writers. Record process/cgroup peak, swaps, temporary
  disk, explicit spill/failure and stable candidate/score/quantity reconciliation.
- [ ] Compatibility: run affected existing TOPP output comparisons and opt-in outputs.
  Preserve supported field semantics and default behavior, including failure modes.
- [ ] Reproducibility: save exact commands, hashes, logs and comparison tolerances.
  Time cold construction and warm reuse end to end at matched search/quality scope.
  Passing tiny fixtures or extraction timing cannot establish greenfield speed gains.

## Risks and completion handoff

- Older binaries versus newer references are not current-source acceptance.
  Resolve build prerequisites and output regressions before claiming validation.
- Baseline mobility/IPF observations do not prove a joint new-engine workflow.
- Hand WP-10 reproducible commands, quality/resource matrices, compatibility results
  and explicit failed/unmeasured gates. Keep the engine opt-in until release gates pass.
