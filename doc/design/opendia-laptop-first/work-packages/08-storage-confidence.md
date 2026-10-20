# WP-08: Localization-capable storage and global confidence

- Status: proposed / not implemented.
- Suggested issue title: `[OpenDIA][WP-08] Compact retained storage and global confidence`.
- Parent epic: OpenDIA laptop-first redesign (issue link pending).
- Issue: not created. This document is a posting draft.
- Repository source: `doc/design/opendia-laptop-first/work-packages/08-storage-confidence.md`.
- Reviewer: unassigned.
- Owner: unassigned. Suggested expertise: statistical calibration, disk-backed I/O.
Repository context: `doc/design/opendia-laptop-first/README.md` and `doc/design/opendia-laptop-first/architecture.md`.

## Goal

Persist sufficient evidence for localization and rescoring independently of the
chosen storage backend, then compute bounded global confidence without combining
independently calibrated tile q-values. Keep identification, localization and
protein-level confidence distinct, including explicit unsupported/unresolved states.

## Dependencies and handoff

- Requires WP-01 contracts.
- Requires WP-02 resource budget.
- Schema and bounded writer work can start in parallel with retrieval/scoring.
- Final integration requires WP-07 family scoring.
- Retrieval-selection assumptions come from WP-06.
- Supplies storage/inference/export contracts to WP-09.

## Scope and non-goals

- Version working evidence, score tables, confidence fields and export mappings.
- Store stable candidate/family/form IDs, shared/identifying fragment attribution,
  ambiguity groups, acquisition/calibration references and rescoring inputs.
- Preserve raw references and query provenance needed for full-evidence revisits.
- Support localization-capable inference without letting a columnar versus other
  backend silently select a scientifically reduced workflow.
- Stream global sorting/reduction and bounded model-fitting passes to disk.
- Define tested populations, competition units and aggregation/error scopes.
- Treat target/decoy selection symmetry as necessary but insufficient validation.
- Do not claim every export is a lossless working store. Declare export capability.
- Do not replace the statistical model solely on an assumed inference bottleneck.

## Deliverables and interface ownership

- Own versioned evidence/score schemas, backend capability declarations and readers.
- Own disk-backed global confidence orchestration and its population metadata.
- Separate identity PEP/q-value, calibrated localization confidence and protein
  confidence, with defined null values and no identity-to-localization substitution.
- Specify stable ordering, boundary deduplication and resumable stage provenance.
- Record exact versus sampled fitting, sample selection, seeds and model versions.
- Define temporary-disk quotas, transactional completion markers and cleanup rules.
- Budget writer queues, encoding scratch, sort buffers, mapped pages and models.
- Reject unsupported requested capabilities before claiming a successful artifact.
- Coordinate existing output adapters with WP-09 and evidence semantics with WP-07.

## Code placement

- Proposed NEW directories, currently absent: `src/openms/include/OpenMS/ANALYSIS/OPENDIA/`
  and `src/openms/source/ANALYSIS/OPENDIA/` own OpenDIA retained schemas,
  retention profiles and inference orchestration. Generic backend machinery belongs
  in `FORMAT`, sorting in `SYSTEM`/`DATASTRUCTURES`, generic math in `MATH/STATISTICS`
  and reusable scoring in `ANALYSIS/OPENSWATH` or appropriate core. No reusable core
  may depend on OPENDIA types/policies. Keep scientific orchestration backend-neutral.

## Output-size contract

- Default compact results preserve ID/localization/quantification ambiguity and run
  provenance with stable IDs/dictionaries and sparse shared evidence, not a duplicate
  full library or dense traces. Lossless defaults must not hide lossy compression.
- Optional rescore evidence retains all declared populations needed by supported
  global rescoring, including rejected targets, decoys and competing alternatives.
  Results-only cannot claim rescoring. Explicit diagnostic traces are selected/capped.
- Profiles change retention, not current-run science. A capability manifest names
  supported future rescore/re-extraction operations and required preserved data/models;
  references alone guarantee neither availability nor reproducibility.
- Separate working inference scratch through global FDR, final artifacts and reusable
  raw/library caches with explicit quotas. Budget peak per-volume totals including
  WAL/journal, indexes, staging and compaction coexistence, not just final file bytes.
- Never drop hypotheses, truncate or alter scientific inputs to meet quotas. Fail
  explicitly before committed final output if mandatory retention cannot fit. Remove
  intermediates only after global inference, validated commit and retention checks.
- Own a concrete schema footprint audit: per-table row fan-out, dictionary/index
  overhead and bytes per meaningful result/evidence, including rejected populations.
  Require declared numerical size goals from WP-00 before evaluation; missing goals
  block size acceptance. Compare OSW at matched capabilities and scientific scope.
- Streaming/Parquet alone is not proof of compactness. Full OSW compatibility export
  is optional, with size preflight and coexistence budget, never the compact default.

## Small PR slices

1. Define versioned schemas and capability tests with synthetic WP-01 records.
   Include ambiguous/coexisting forms and backend-independent inference inputs.
2. Add bounded append/read, external sort/reduction and interruption/disk-full
   handling. Verify round trips before introducing confidence calculations.
3. Implement declared global competition and bounded fitting/calibration passes,
   distinguishing exact and sampled execution and separate confidence scopes.
4. Integrate real WP-07 records, localization inference and existing output adapters.
   Add cross-backend and public CLI tests without changing default workflows.

## Requirement-to-acceptance checks

- [ ] Review includes and build dependencies: retained schema/profile/inference policy
  stays in OPENDIA; generic backend/sorting code has no OPENDIA types or policies.
- [ ] Audit scaled-library/family fixtures against WP-00 numerical size goals and
  matched-capability OSW, reporting row fan-out and bytes/result/evidence by profile.
- [ ] Reopen committed profiles, verify manifest/ambiguity/provenance and full rescore
  populations; reject unsupported operations and test missing referenced data/models.
- [ ] Inject quota/disk-full/interruption during WAL/index/staging/compaction and OSW
  export; verify peak per-volume accounting, commit safety and retention-aware cleanup.
- [ ] Storage capability: round-trip shared and identifying evidence, IDs, ambiguity,
  provenance and missing-value semantics through each advertised working backend.
  Rescoring must retain localization capability and equivalent scientific results.
- [ ] Global FDR: compare disk-backed output to an in-memory small-case global oracle.
  Changing tile partition/order must not change the tested population or q-values
  beyond predeclared numeric tolerance. Never concatenate tile-level estimates.
- [ ] Selection validity: evaluate held-out/entrapment and controlled absent hypotheses
  across abundance, modification type and ambiguity, including retrieval pruning.
  Equivalent target/decoy filtering alone does not pass the calibration gate.
- [ ] Localization: known-site mixtures and unresolved forms validate localization
  error separately from identification. Protein aggregation has its own scope.
- [ ] Fitting: exact/sampled modes report their choice and reproduce seeded results.
  Compare sampled calibration to exact fitting on tractable held-out fixtures.
- [ ] Bounded storage: grow rows/absent hypotheses and throttle consumers at 2/4 GiB.
  Measure process peak and disk usage, and test quota, interruption and disk-full
  failures without partial outputs being mistaken for completed calibrated results.
- [ ] Public acceptance after WP-09: actual opt-in OpenDIA CLI writes stored evidence
  from raw mobility input, then supported public inference/export entry points
  consume it for localization. Compare supported backends and final result fields.

## Risks and completion handoff

- Grouping, recovery and competition change the tested population. Freeze and
  validate those semantics before treating nominal q-values as calibrated errors.
- Existing Parquet/IPF rejection is a baseline limitation, not an acceptable silent
  downgrade for the new engine or proof that a new schema already works.
- Hand off schema versions, capability matrix, calibration reports, disk estimates
  and restart/error behavior to WP-09. Report unresolved statistical limitations.
- These are implementation gates, not measured greenfield performance claims.
