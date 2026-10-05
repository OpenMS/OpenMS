# OpenMS identification data: final design

**Status:** simplified native format and owning model, reconstructed on 5 October 2026. The implementation uses schema version 1; it has not been released. Validation results and migration limits are listed below.

## 1. Decision

Use an ordinary directory containing a small JSON manifest and typed Parquet tables. Keep analysis runs logically separate while sharing physical tables across compatible runs. Keep scores in the match table. Write a new, self-contained dataset after filtering or modification.

There are no run revision tokens, scientific snapshots, file catalogues, update overlays, persistent lookup indexes or partial-replacement protocol. Assign a public schema version when the format is ready; discarded, unreleased prototypes do not need compatibility readers.

The format supports ordinary editing and sequential processing of datasets larger than RAM. It does not promise that a billion complete C++ match objects, or arbitrary global inference, fit in memory.

## 2. Data model and run contract

An **analysis run** groups identifications with compatible configuration and score definitions. It can contain observations from several physical MS files. A **source** records the originating file or an explicitly unknown source within a run; it can retain an ordered primary-file list. Source identity is never inferred from a basename.

Each run defines:

- Its identifier, molecule kind and search/processing metadata.
- Its sources, including exact paths and optional primary-file lists.
- Immutable score definitions and an optional primary score.
- Shared metadata descriptors: name, value type and unit information.

A query represents an observation with zero or more candidate matches. It has a generic native reference such as a spectrum identifier, optional RT and observed m/z. There is no mandatory spectrum-versus-feature type flag. Feature associations are separate annotations.

A match owns its molecular representation, ion description, scores, original parent evidence and optional annotations. Peptides, oligonucleotides and compounds use the same structure. Molecule kind is a run constant; representation is a string with an explicit encoding, such as AASequence notation, NASequence notation, SMILES or a database identifier. Ordinary I/O does not parse or canonicalize it.

## 3. Directory layout

The manifest identifies every run/result and its table slices: relative physical path, starting row, row count, and numeric partition ID. The partition maps to the run/result descriptor and its stable UUID; it is a storage coordinate, not a new scientific identity. Human identifiers and source paths are never used directly as filenames.

| Path | Contents |
| --- | --- |
| `manifest.json` | Format identifier/schema version, run descriptors, sources, score definitions, metadata descriptors, processing metadata and table paths |
| `queries-0.parquet` | Queries for all runs |
| `matches-0.parquet` | Matches and all their score columns |
| `parents-0.parquet` | Optional original parent catalogue: qualified identity, sequence, description and metadata |
| `inputs-0.parquet` | Ordered contributing runs and exact input score definitions |
| `input_members-0.parquet` | Ordered candidate memberships |
| `proteins-0.parquet` | Inferred protein hits |
| `groups-0.parquet` | Protein group values and ordering |
| `group_members-0.parquet` | Ordered members of each group |
| `assignments-0.parquet` | Ordered match assignments with a typed parent list, including empty lists |

Use one physical file per compatible table schema, with multiple Parquet row groups. Small runs and inference results share row groups. Query, match and parent rows carry `run_id`; inference rows carry `inference_id`. Different complete supplementary score layouts use separate match files (`matches-1.parquet`, etc.), preserving dense numeric columns and their definitions. Files do not multiply with the number of compatible runs. Size-based file sharding remains outside this initial implementation.

Keep JSON for configuration whose size follows the number of runs, sources, score definitions and metadata descriptors. Per-record values, parent catalogues, protein results, memberships and assignments belong in typed tables. The manifest contains no arrays of PSM IDs.

The manifest stores small inference-result descriptors: result ID, display identifier, algorithm/parameters, output score definitions, metadata descriptors and table paths. Detailed per-protein/group metadata is stored with its table rows.

Every run declares query and match slices, even when they contain zero rows. An omitted parent table means no parent catalogue was supplied; a declared zero-row table means a supplied catalogue was empty. Every inference result declares slices in the six shared tables above; unused collections have zero-row slices. An identification-only dataset has an empty result list and no inference tables. Missing declared files are errors.

## 4. Core tables

### Queries

| Field | Meaning |
| --- | --- |
| `query_id` | Nonzero uint64 identifier within the run |
| `source_id` | uint32 reference to a source descriptor |
| `data_id` | Exact original observation reference |
| `rt`, `mz` | Nullable float64; RT in seconds, observed m/z in Th |
| `selected_match_id` | Optional explicitly selected candidate |
| `metadata` | Typed metadata values |

### Matches

| Field | Meaning |
| --- | --- |
| `match_id` | Nonzero uint64 identifier within the run |
| `query_id` | Owning query |
| `representation`, `encoding` | Exact molecular string and its interpretation |
| `charge`, `calculated_mz` | int32 charge and nullable float64 m/z in Th |
| `target_decoy` | Explicit unknown, target, decoy or both state |
| `name`, `formula`, `identifiers` | Molecular name, optional formula and qualified identifiers |
| `adduct` | Nullable struct containing name, formula, net charge and multiplier |
| `parent_evidence` | Ordered sequence-to-parent evidence with qualified identities, optional coordinates and flanking residues |
| `peak_annotations` | Ordered fragment annotations |
| `metadata` | Typed metadata values |
| `score_<id>` | One nullable float64 column per run score definition |

Queries and candidates retain their order. Match rows follow query order, with a query's candidates contiguous. A reader advances through the query stream and consumes consecutive match rows whose query ID equals the current query, retaining one row of lookahead. A query with no matching row is empty. Any unmatched match rows after the query stream ends are invalid. Validate the selected candidate while traversing its query. This needs neither a stored candidate count nor a run-wide query-to-match hash map. Queries with no candidates remain representable. A very large candidate set can span row groups and processing batches.

The selected match is independent of the primary score and any calculated ranking. Evidence positions are zero-based and inclusive, with independently nullable endpoints. Adduct net charge must agree with ion charge; multiplier must be positive. Known parent identities include their database namespace. Compounds do not acquire fabricated sequence-to-parent evidence.

## 5. Scores and repeated information

A score definition includes name, optional CV term, direction, statistical scope, producing software/parameters and calibration provenance when available. Equal names alone do not establish equal meanings. The inference algorithm is responsible for checking whether scores from different runs can be pooled.

A score has a stable run-local uint32 ID. The complete definition set is known before opening the match writer, so every row group has the same score schema. Every `score_<id>` column has exactly one meaning throughout its run. Missing scores use null; zero remains a value. New-model scores are finite. Selecting a primary score requires a value on every match, including subsequently appended matches.

Changing the primary score uses the dataset-wide checked selection API. Rescoring with a different meaning or provenance creates a new definition and column. Inference retains the definition used for its calculation.

Keep score columns in the shared match tables. Parquet column projection supports score-only reads without decoding molecule and annotation columns. Separate score files would introduce alignment rules and are not required for this operation. Rewriting the dataset after rescoring is an accepted first-version tradeoff.

Source paths appear once per source descriptor, not on every query or PSM. Repetition of a path across a few run descriptors is acceptable. Optional adducts remain owned nullable structs; a global adduct registry is unnecessary.

Retain the prototype's typed metadata codec and shared descriptors. Preserve differences between an absent key, EMPTY_VALUE, an empty string and an empty typed list. Preserve value types, units, list order and supported floating-point bit patterns. Do not turn typed values into JSON strings.

The descriptor dictionary and manifest remain resident while scanning. Their memory cost must be reported separately from batch buffers; datasets with unusually many distinct descriptors do not have a constant metadata footprint. The initial native text contract is valid UTF-8, without normalization or trimming. Writers reject unsupported byte strings rather than silently converting them. Exact preservation of arbitrary non-UTF-8 legacy strings is outside this initial contract. Native exactness is measured against the supported owning model; legacy adapters report or reject unrepresentable fields.

## 6. Stable identities and filtering

Each analysis run has an opaque UUID plus a human-readable identifier. Persistent query/match identity is `(run_uuid, record_kind, record_id)`. Query and match IDs have separate uint64 namespaces. Surviving IDs are preserved when writing a subset; row positions can change.

The manifest retains the next allocatable query and match IDs for each run, including across filtering. These are ID allocation counters, not revisions: score or metadata edits do not advance them. Deleted IDs are never reused. For a present run, the counter must exceed IDs in live rows and preserved provenance; importing an externally supplied ID reserves it as well. An absent run UUID remains unresolved and cannot be reused for an independent run. Exhaustion is an error.

A copied subset keeps its run identity so existing provenance still identifies its records. A new independent analysis gets a new run UUID. Merging overlapping or divergent copies requires an explicit conflict decision; the importer must reject ambiguity or remap identities and all associated links. It must not silently treat equal IDs as equal payloads.

Filtering processes a run sequentially and writes the retained queries and matches to a new directory. It clears a selection if the selected match was removed. Filtering removes empty queries by default; `keep_empty_queries=true` retains them. Saving without filtering always preserves existing empty queries. No renumbering of surviving matches is needed.

Editing metadata, scores or candidate payloads uses the owning API's validity rules. Saving writes the resulting values. Native file IDs and process-local handles need not be identical: streaming views use file IDs directly, while an owning adapter may allocate fresh runtime handles and return mappings.

## 7. Pooled inference

Inference results belong to the dataset, independently of individual runs. One result can combine any number of analysis runs and stores its protein/group output once.

The six inference files use result-local keys and explicit ordinals. Reuse the supported OpenMS protein field definitions and typed metadata codecs. Large membership and group-member collections are rows in child tables, never one dataset-sized list cell. Each assignment is one row with its ordered typed parent list; this list has the same individual-payload memory consideration as a match's original parent evidence.

| Table | Required content |
| --- | --- |
| `inputs` | Input ID/order, run UUID, optional exact score definition, selection provenance, membership-known flag, member count when known |
| `input_members` | Input ID, member ordinal, match ID; preserve order and duplicates |
| `proteins` | Result-local protein key/order, qualified identity, original alias, supported ProteinHit values and typed metadata |
| `groups` | Group key/order, general or indistinguishable kind, original score with declared meaning, member count |
| `group_members` | Group key, member ordinal, original alias and optional qualified identity; do not require an invented protein hit |
| `assignments` | Assignment key/order, run UUID, match ID, optional input ID and ordered list of qualified parents |

Child rows are grouped by their input or group header and follow its recorded order/count. They can be consumed sequentially without an inference-wide hash join.

Membership means the candidates considered by the calculation; it is distinct from winning candidates and protein assignments. A known empty membership differs from unknown membership. An assignment row with an empty parent list differs from an absent assignment. Missing membership rows never mean all current matches. If an imported result supplies only run-level provenance, preserve that limitation rather than inventing exact membership.

Dataset-level transformations require the workflow to choose the treatment of affected inference results:

1. **Preserve:** keep the original inference output, scores, memberships and assignments.
2. **Discard:** remove the inference output and write identification-only data.
3. **Recompute:** discard the old result and explicitly rerun the chosen algorithm.

Low-level run edits preserve inference and do not walk inference objects. Dataset-level transformations receive a policy argument; no interactive confirmation is needed. The serializer performs structural validation and writes the selected contents; it does not judge scientific freshness or run inference automatically.

Preserved inference may mention removed candidates or absent input runs. Those IDs describe the original calculation and are allowed to remain unresolved. They must remain distinct from live IDs during loading. A surviving ID may also refer to a subsequently edited hypothesis, so ID existence does not establish that an old assignment applies to its current contents.

No old molecular payloads are retained. Reproducing an earlier calculation requires independently saved inputs.

Filtering inferred proteins keeps surviving protein values, removes assignments to removed parents, and retains an explicit empty parent list when an assignment's last parent disappears. By default it drops every group containing a removed member. Recalculation is an explicit algorithm call. The original parent catalogue and search evidence remain unchanged.

## 8. Reading and writing

Provide two access paths with the same semantics:

- **Owning access:** load a run or complete dataset for conventional OpenMS algorithms. Object allocation and optional lookup maps consume memory proportional to the loaded data.
- **Streaming access:** select runs and columns, decode batches, and write results incrementally. Sequential scans and filters require no match-sized lookup map.

Use configurable Arrow batch and Parquet row-group sizes with both byte targets and row-count limits. A filtering writer can stream a query's matches first, then emit its query row with the surviving selection; it does not buffer all candidates for a large query. Start with the current compression settings and tune through Release benchmarks. A single unusually large payload remains an explicit memory consideration.

Selecting one of 1,000 runs opens only that run's requested tables after reading the manifest. Billion-PSM processing relies on streaming or loading manageable runs. Global inference may still need substantial memory or a dedicated external-memory algorithm; file organization alone does not solve that.

Write a fresh output directory through a temporary sibling directory. Close and validate files, then publish the completed directory. Reject an existing destination in the first implementation. This avoids introducing concurrent updates, table replacement, transaction catalogues and retained file generations. Atomic publication and power-loss durability are separate platform-dependent guarantees.

## 9. Validation and failure behavior

Writers enforce the run contract, ID allocation rules and supported value types. Readers validate schemas, source/score descriptors, row counts, typed payloads and sequential ownership. Selected matches must occur within their query. Inference memberships and assignments may retain unresolved provenance IDs; live feature links must resolve or be explicitly external.

ID uniqueness is a separate global check. A normal sequential scan needs no record-sized map. Strict validation of untrusted/imported IDs may use a set for a manageable run or external sorting for a large run; it must preserve the stored scientific order. Random lookup likewise builds an explicit optional index in memory. Neither cost is hidden in the streaming guarantee.

Owning load validates into a temporary object and preserves the destination on failure. A streaming reader can report a late error after earlier batches were delivered. A failed writer leaves no published output dataset.

Custom modification definitions are retained with run configuration. Raw scans preserve molecular strings without registering global chemistry. Explicit chemical materialization resolves definitions in the run context and rejects conflicting unsupported definitions rather than reusing an unrelated global entry.

## 10. OpenMS workflow implications

| Workflow | Intended use |
| --- | --- |
| IDScoreSwitcher | Resolve a run score definition and select it as primary; validate coverage |
| Protein inference | Read selected scores and evidence across runs; write one pooled result with explicit inputs |
| NASE / ProSE | Use homogeneous RNA or peptide runs with the appropriate molecular encoding and score contract |
| ProteomicsLFQ | Consume identifications and retain explicit associations to measured features |
| IsobaricWorkflow | Retain quantitative channels, feature associations and pooled inference as distinct information |

Feature and consensus persistence can reuse these identification and inference tables through optional association tables. Full quantification/map schemas are a later implementation step, not extra machinery required by the base identification format. Removing a PSM does not remove a measured feature; live associations to deleted PSMs must be pruned or the operation rejected. These live links differ from retained inference provenance.

## 11. Implementation and migration

The implementation consists of:

- `IdentificationData`: owned runs, source blocks, queries, candidates and independent inference results. Numeric scores are stored densely, with missing values exposed as optionals. Runtime lookup maps are lazy; call `prepareLookupIndexes()` before sharing a run for parallel read-only lookups.
- `IdentificationDataFile`: manifest plus typed Parquet tables, stable persisted identities, multi-row-group writing, projected scans, run loading and streaming filtering into a fresh dataset.
- `IdentificationDataAdapter`: legacy peptide/protein conversion and feature/consensus association adapters. Strict conversion rejects information the destination model cannot express; explicit permissive conversion reports losses.
- `IdentificationDataInference`: a pooled BasicProteinInference bridge with explicit probability score selection, qualified parent identities, membership and assignment provenance, plus protein filtering.
- Python bindings for the owning model, adapters, inference, native I/O and FileHandler access.

The former reference-based implementation is named `LegacyIdentificationData` during consumer migration. Existing OMS, feature internals, NASE and other algorithms still using that representation are updated mechanically to the explicit legacy name. They are not implicitly converted into the new model. Feature and consensus adapters enable deliberate migration while retaining quantitative measurements and checking live links. Migrating every workflow and persisting complete quantitative maps are separate follow-up tasks.

`FileHandler` recognizes the native manifest and provides owning load/store overloads. Legacy idXML input can be imported; legacy export uses strict adapters. The established PSM Parquet bundle remains a distinct supported input and is not confused with the native format. No reader for discarded experimental versions is provided.

The benchmark in `tools/benchmarks/identification_data` measures owning generation/load, native writing, score-only scans, full scans and streaming filtering in separate Release processes. Run count and match count are configurable. Synthetic scaling and real idXML measurements must be reported separately. Descriptor residency, Arrow buffers and individual large payloads are additional to callback batch targets; byte targets are not a hard process-memory limit.

See the benchmark README for reproducible commands and the accompanying validation report for measured results.

## Checked common primary PSM score

A dataset has one common primary PSM score definition across participating runs.
The complete definition must match, including direction, scope and provenance;
matching display names alone is insufficient. Local score IDs and supplementary
score columns may differ. Protein and inference scores remain independent.
Empty runs with no selected primary are ignored during incremental construction.
Populated runs must either all select the same definition or all be unscored.

`getPrimaryScoreDefinition()` checks and returns this contract.
`setPrimaryScore(definition)` checks score availability and candidate coverage in
every configured run before changing any selection. Failure changes no selections.
`addRun`, `replaceRun`, legacy import, native read/write, and inference reject
incompatible datasets. Native descriptor checks also cover streaming scans.
Mutable C++ run access supports incremental edits; call `validate()` after such
edits. The contract is checked at these boundaries, not on every low-level setter.
Normalize mixed legacy scores before import; splitting them into runs no longer
makes them a valid single dataset. Supplementary schemas may differ; the shared-file writer groups compatible complete score layouts.

## Shared-table I/O implementation

The writer fills Arrow builders directly from references to the owning query/match
values and existing numeric score vectors. It creates no temporary owning payloads
or optional-score vectors per PSM. Logical table slices close without flushing;
physical builders flush at the configured row/byte targets and once at final close.

Readers select intersecting row groups using manifest row ranges, trim them to the
requested slice, and validate the partition column. Projected decoding caches one
row group per open physical table, with up to 16 idle/active file entries retained
when possible (active readers are never evicted). Sequential runs sharing a row group
reuse its decoded arrays. Callback batches still obey their row/estimated byte
limits; cached decoded row groups use additional memory bounded by physical row
group size, not callback batch size. External oversized row groups can exceed those
writer targets. Projection excludes unused payload columns.

Filtering rewrites shared query/match tables and their row ranges. Retained parent
and inference tables are validated by slice and copied once per physical file.
Empty versus absent catalogues, local IDs, ordering and inference provenance retain
their previous semantics. Malformed overlapping ranges or mismatched partitions
are rejected. No reader for the discarded unreleased per-run-file layout is kept.
