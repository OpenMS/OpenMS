# WP-08: Localization-capable storage and global confidence

- Status: proposed / not implemented.
- Reviewer: unassigned.
- Owner: unassigned. Suggested expertise: statistical calibration, disk-backed I/O.
- Context: [master plan](../README.md), [architecture](../architecture.md).

## Goal

Persist sufficient evidence for localization and rescoring independently of the
chosen storage backend, then compute bounded global confidence without combining
independently calibrated tile q-values. Keep identification, localization and
protein-level confidence distinct, including explicit unsupported/unresolved states.

## Dependencies and handoff

- Requires [WP-01 contracts](01-contracts.md).
- Requires [WP-02 resource budget](02-resource-budget.md).
- Schema and bounded writer work can start in parallel with retrieval/scoring.
- Final integration requires [WP-07 family scoring](07-family-scoring.md).
- Retrieval-selection assumptions come from [WP-06](06-candidate-retrieval.md).
- Supplies storage/inference/export contracts to [WP-09](09-cli-integration.md).

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

- Storage capability: round-trip shared and identifying evidence, IDs, ambiguity,
  provenance and missing-value semantics through each advertised working backend.
  Rescoring must retain localization capability and equivalent scientific results.
- Global FDR: compare disk-backed output to an in-memory small-case global oracle.
  Changing tile partition/order must not change the tested population or q-values
  beyond predeclared numeric tolerance. Never concatenate tile-level estimates.
- Selection validity: evaluate held-out/entrapment and controlled absent hypotheses
  across abundance, modification type and ambiguity, including retrieval pruning.
  Equivalent target/decoy filtering alone does not pass the calibration gate.
- Localization: known-site mixtures and unresolved forms validate localization
  error separately from identification. Protein aggregation has its own scope.
- Fitting: exact/sampled modes report their choice and reproduce seeded results.
  Compare sampled calibration to exact fitting on tractable held-out fixtures.
- Bounded storage: grow rows/absent hypotheses and throttle consumers at 2/4 GiB.
  Measure process peak and disk usage, and test quota, interruption and disk-full
  failures without partial outputs being mistaken for completed calibrated results.
- Public acceptance after WP-09: actual opt-in OpenDIA CLI writes stored evidence
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
