# OpenMS identification data: final design

**Status:** simplified native format and owning model, reconstructed on 5 October 2026. The implementation uses schema version 1; it has not been released. Validation results and migration limits are listed below.

## 1. Decision

Use an ordinary directory containing a small JSON manifest and typed Parquet tables. Keep analysis runs logically separate while sharing physical tables across runs with one checked score schema. Keep scores in the match table. Write a new, self-contained dataset after filtering or modification.

There are no run revision tokens, scientific snapshots, file catalogues, update overlays, persistent lookup indexes or partial-replacement protocol. Assign a public schema version when the format is ready; discarded, unreleased prototypes do not need compatibility readers.

The format supports ordinary editing and sequential processing of datasets larger than RAM. It does not promise that a billion complete C++ match objects, or arbitrary global inference, fit in memory.

## 2. Data model and run contract

An **analysis run** groups identifications with compatible configuration and score definitions. It can contain observations from several physical MS files. A **source** records the originating file or an explicitly unknown source within a run; it can retain an ordered primary-file list. Source identity is never inferred from a basename.

Each run defines:

- Its identifier, molecule kind and search/processing metadata.
- Its sources, including exact paths and optional primary-file lists.
- The shared ordered score definitions and primary selection; empty unconfigured runs may omit both.
- Shared metadata descriptors: name, value type and unit information.

A query represents an observation with zero or more candidate matches. It has a generic native reference such as a spectrum identifier, optional RT and observed m/z. There is no mandatory spectrum-versus-feature type flag. Feature associations are separate annotations.

A match owns its molecular representation, ion description, scores, original parent evidence and optional annotations. Peptides, oligonucleotides and compounds use the same structure. Molecule kind is a run constant; representation is a string with an explicit encoding, such as AASequence notation, NASequence notation, SMILES or a database identifier. Ordinary I/O does not parse or canonicalize it.

## 3. Directory layout

The manifest identifies every run/result and its table slices: relative physical path, starting row, row count, and numeric partition ID. The partition maps to the run/result descriptor and its stable UUID; it is a storage coordinate, not a new scientific identity. Human identifiers and source paths are never used directly as filenames.

| Path | Contents |
| --- | --- |
| `manifest.json` | Format identifier/schema version, run descriptors, sources, score definitions, metadata descriptors, processing metadata and table paths |
| `queries.parquet` | Queries for all runs |
| `matches.parquet` | Matches and all their score columns |
| `parents.parquet` | Optional original parent catalogue: qualified identity, sequence, description and metadata |
| `inputs.parquet` | Ordered contributing runs and exact input score definitions |
| `proteins.parquet` | Inferred protein hits |
| `groups.parquet` | Protein group values, ordering and an inline ordered member list |

The layout above has at most seven files including the manifest; empty optional tables need no physical file. Use exactly one physical file per table, with multiple Parquet row groups. Small runs and inference results share row groups. Query, match and parent rows carry `run_id`; inference rows carry `inference_id`. Every configured run uses the same ordered score definitions and primary column. Supplementary values may be null. Different score layouts are rejected instead of creating additional files. Filenames have no numeric suffix. Size-based file sharding remains outside this initial implementation.

Keep JSON for configuration whose size follows the number of runs, sources, score definitions and metadata descriptors. Per-record values, parent catalogues, protein results and group members belong in typed tables. The manifest contains no arrays of PSM IDs.

The manifest stores small inference-result descriptors: result ID, display identifier, algorithm/parameters, output score definitions, metadata descriptors and table paths. Detailed per-protein/group metadata is stored with its table rows.

Every run declares query and match slices, even when they contain zero rows. An omitted parent table means no parent catalogue was supplied; a declared zero-row table means a supplied catalogue was empty. Every inference result declares slices in the three shared inference tables above; unused collections have zero-row slices. An identification-only dataset has an empty result list and no inference tables. Missing declared files are errors.

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
| `score_<id>` | One float64 column per shared score definition; supplementary values are nullable |

Queries and candidates retain their order. Match rows follow query order, with a query's candidates contiguous. A reader advances through the query stream and consumes consecutive match rows whose query ID equals the current query, retaining one row of lookahead. A query with no matching row is empty. Any unmatched match rows after the query stream ends are invalid. Validate the selected candidate while traversing its query. This needs neither a stored candidate count nor a run-wide query-to-match hash map. Queries with no candidates remain representable. A very large candidate set can span row groups and processing batches.

The selected match is independent of the primary score and any calculated ranking. Evidence positions are zero-based and inclusive, with independently nullable endpoints. Adduct net charge must agree with ion charge; multiplier must be positive. Known parent identities include their database namespace. Compounds do not acquire fabricated sequence-to-parent evidence.

## 5. Scores and repeated information

A score definition includes name, optional CV term, direction, statistical scope, producing software/parameters and calibration provenance when available. Equal names alone do not establish equal meanings. The inference algorithm is responsible for checking whether scores from different runs can be pooled.

A score has a stable uint32 column index shared by every configured run; C++ ScoreId handles remain bound to their owning run. The complete ordered definition set is checked before opening the match writer, so every row group has the same score schema. Every `score_<id>` column has exactly one meaning throughout the dataset. Missing supplementary scores use null; zero remains a value. New-model scores are finite. Selecting a primary score requires a value on every match, including subsequently appended matches.

Changing the primary score uses the dataset-wide checked selection API. Rescoring with a different meaning or provenance creates a new definition and column. Inference retains the definition used for its calculation.

Keep score columns in the shared match tables. Parquet column projection supports score-only reads without decoding molecule and annotation columns. Separate score files would introduce alignment rules and are not required for this operation. Rewriting the dataset after rescoring is an accepted first-version tradeoff.

Source paths appear once per source descriptor, not on every query or PSM. Repetition of a path across a few run descriptors is acceptable. Optional adducts remain owned nullable structs; a global adduct registry is unnecessary.

Retain the prototype's typed metadata codec and shared descriptors. Preserve differences between an absent key, EMPTY_VALUE, an empty string and an empty typed list. Preserve value types, units, list order and supported floating-point bit patterns. Do not turn typed values into JSON strings.

The descriptor dictionary and manifest remain resident while scanning. Their memory cost must be reported separately from batch buffers; datasets with unusually many distinct descriptors do not have a constant metadata footprint. The initial native text contract is valid UTF-8, without normalization or trimming. Writers reject unsupported byte strings rather than silently converting them. Exact preservation of arbitrary non-UTF-8 legacy strings is outside this initial contract. Native exactness is measured against the supported owning model; legacy adapters report or reject unrepresentable fields.

## 6. Stable identities and filtering

Each analysis run has an opaque UUID plus a human-readable identifier. Persistent query/match identity is `(run_uuid, record_kind, record_id)`. Query and match IDs have separate uint64 namespaces. Surviving IDs are preserved when writing a subset; row positions can change.

The manifest retains the next allocatable query and match IDs for each run, including across filtering. These are ID allocation counters, not revisions: score or metadata edits do not advance them. Deleted IDs are never reused. For a present run, the counter must exceed IDs in live rows; importing an externally supplied ID reserves it as well. An absent run UUID remains unresolved and cannot be reused for an independent run. Exhaustion is an error.

A copied subset keeps its run identity so existing provenance still identifies its records. A new independent analysis gets a new run UUID. Merging overlapping or divergent copies requires an explicit conflict decision; the importer must reject ambiguity or remap identities and all associated links. It must not silently treat equal IDs as equal payloads.

Filtering processes a run sequentially and writes the retained queries and matches to a new directory. It clears a selection if the selected match was removed. Filtering removes empty queries by default; `keep_empty_queries=true` retains them. Saving without filtering always preserves existing empty queries. No renumbering of surviving matches is needed.

Editing metadata, scores or candidate payloads uses the owning API's validity rules. Saving writes the resulting values. Native file IDs and process-local handles need not be identical: streaming views use file IDs directly, while an owning adapter may allocate fresh runtime handles and return mappings.

Retain IDs across edits and reacquire references to match records after filtering, replacement or transformation. These operations may replace query-owned match storage even when record counts do not change. Edits stage potentially throwing work before committing; platforms with nonthrowing payload moves retain the in-place fast path. Run processing metadata is owned independently and replaced through a nonthrowing pointer swap.

## 7. Pooled inference

Inference results belong to the dataset, independently of individual runs. One result can combine any number of analysis runs and stores its protein/group output once.

The three inference files use result-local keys and explicit ordinals. Reuse the supported OpenMS protein field definitions and typed metadata codecs. Input provenance is recorded per contributing run; it does not contain a separate list of candidate IDs. Each protein group stores its ordered members directly in an Arrow list of structs (original alias and optional qualified identity). List order preserves duplicates and an empty list preserves an empty group. Group membership needs no separate file, join key, ordinal or member-count column. A complete group row is subject to `max_record_bytes`.

| Table | Required content |
| --- | --- |
| `inputs` | Input ID/order, run UUID and identifier, optional exact score definition, selection description |
| `proteins` | Result-local protein key/order, qualified identity, original alias, supported ProteinHit values and typed metadata |
| `groups` | Group key/order, general or indistinguishable kind, original score, typed arrays, and `members: list<struct<alias, identity?>>`; group-only aliases do not require an invented protein hit |

Input records describe which analysis runs and scores contributed, plus a human-readable selection description. They do not assert that all current matches were used, and the exact candidate set cannot be reconstructed from these records. Reproducing the calculation requires independently saved inputs and selection parameters.

The base model and format do not store inferred PSM-to-protein assignments. Original search evidence remains on each match. A future inference model may represent PSM-to-peptide-to-protein relationships in a separate format; its graph structure and persistence are outside this design.

If a workflow needs membership annotations for retained matches, it can use optional match metadata. Such annotations disappear with deleted matches and do not form a retained inference graph.

Dataset-level transformations require the workflow to choose the treatment of affected inference results:

1. **Preserve:** keep the original inference output, scores and run-level input provenance.
2. **Discard:** remove the inference output and write identification-only data.
3. **Recompute:** discard the old result and explicitly rerun the chosen algorithm.

Low-level run edits preserve inference and do not walk inference objects. Dataset-level transformations receive a policy argument; no interactive confirmation is needed. The serializer performs structural validation and writes the selected contents; it does not judge scientific freshness or run inference automatically.

Preserved inference may refer to absent input runs. Their UUIDs describe the original calculation and are allowed to remain unresolved. Run identity alone does not establish that an inference result was calculated from the currently stored matches.

No old molecular payloads are retained. Reproducing an earlier calculation requires independently saved inputs.

Filtering inferred proteins keeps surviving protein values. By default it drops every group containing a removed member. Recalculation is an explicit algorithm call. The original parent catalogue and search evidence remain unchanged.

## 8. Reading and writing

Provide two access paths with the same semantics:

- **Owning access:** load a run or complete dataset for conventional OpenMS algorithms. Object allocation and optional lookup maps consume memory proportional to the loaded data.
- **Streaming access:** select runs and columns, decode batches, and write results incrementally. Sequential scans and filters require no match-sized lookup map.

Use configurable Arrow batch and Parquet row-group sizes with both byte targets and row-count limits. A filtering writer can stream a query's matches first, then emit its query row with the surviving selection; it does not buffer all candidates for a large query. Start with the current compression settings and tune through Release benchmarks. A single unusually large payload remains an explicit memory consideration.

`IdentificationDataFile::Options::threads` is a positive per-operation CPU worker
limit (default 1), also exposed as `Options.threads` in Python. Higher counts use
a private Arrow pool for column decoding and buffered row-group encoding. They do
not alter Arrow's global pool, parallelize model mutations or create more files.
Decoding retains one projected row group per cached physical table with no
row-group readahead; parallel encoding can require additional per-column buffers.
The caller waits for each batch before reusing its storage, and failures propagate
through the existing transactional load/store boundaries. Work on different
operations uses separate pools.

Owning loading transfers decoded payloads by move and validates each completed run
once. Validation reads dense score storage directly and checks uniqueness with
compact ID vectors, using a linear path for sorted IDs and sorting copies otherwise.
Scientific record ordering remains unchanged. Column indices are bound once per
batch and metadata descriptors cache runtime registry IDs; those runtime IDs are
never persisted in the file format.

Selecting one of 1,000 runs opens only that run's requested tables after reading the manifest. Billion-PSM processing relies on streaming or loading manageable runs. Global inference may still need substantial memory or a dedicated external-memory algorithm; file organization alone does not solve that.

Write a fresh output directory through a temporary sibling directory. Close and validate files, then publish the completed directory. Reject an existing destination in the first implementation. This avoids introducing concurrent updates, table replacement, transaction catalogues and retained file generations. Atomic publication and power-loss durability are separate platform-dependent guarantees.

## 9. Validation and failure behavior

Writers enforce the run contract, ID allocation rules and supported value types. Readers validate schemas, source/score descriptors, row counts, typed payloads and sequential ownership. Selected matches must occur within their query. Inference inputs may refer to absent runs; live feature links must resolve or be explicitly external.

ID uniqueness is a separate global check. A normal sequential scan needs no record-sized map. Strict validation of untrusted/imported IDs may use a set for a manageable run or external sorting for a large run; it must preserve the stored scientific order. Random lookup likewise builds an explicit optional index in memory. Neither cost is hidden in the streaming guarantee.

Owning load validates into a temporary object and preserves the destination on failure. A streaming reader can report a late error after earlier batches were delivered. A failed writer leaves no published output dataset.

Custom modification definitions are retained with run configuration. Raw scans preserve molecular strings without registering global chemistry. Explicit chemical materialization resolves definitions in the run context and rejects conflicting unsupported definitions rather than reusing an unrelated global entry.

## 10. OpenMS workflow implications

| Workflow | Intended use |
| --- | --- |
| IDScoreSwitcher | Select a shared score as primary across runs; validate coverage atomically |
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
- `IdentificationDataInference`: a pooled BasicProteinInference bridge with explicit probability score selection, qualified parent identities and run-level input provenance, plus protein filtering.
- Python bindings for the owning model, adapters, inference, native I/O and FileHandler access.

The reference-based implementation and its internal graph classes have been removed. FeatureMap and ConsensusMap own IdentificationData directly. Feature annotations contain qualified match/query IDs and an optional molecular identity value. Copying maps preserves those IDs; it does not translate references. Dataset-aware annotation statistics resolve the owning match values. The no-argument annotation-state API still handles legacy peptide annotations and one native match; comparing several native matches requires the dataset argument.

Filtering, pooled FDR, retention-time alignment, accurate-mass search, NASE and conversion use the owning APIs. Independent runs with repeated display names receive numeric suffixes during merging; their UUIDs and feature links stay unchanged. Conflicting values for the same UUID fail atomically. Protein inference retains explicit run-level provenance across filtered data.

OMS schema 6 keeps quantitative measurements in SQLite and embeds the materialized native identification files as bounded binary chunks. It is a single .oms file and shares the native serializer rather than maintaining another PSM schema. Compatible released OMS schemas are read directly into owning records using temporary import maps. This compatibility path retains historical processing and score applications as typed metadata; the materialized score is its latest application. Opposite or disjoint PSM score layouts without a complete common primary column, and groupings with multiple score types that cannot fit one explicit inference score contract, are rejected. JSON export is a diagnostic view of SQLite, not the native identification representation.

RNA digestion returns owning sequence/evidence candidates. NASE can persist a scoreless catalog run whose observations identify catalog entries. Such runs explicitly set `identification:catalog=true` in processing metadata and cannot declare PSM score columns; ordinary search runs still require a common primary score. RNA/compound idXML export retains the historical label/molecule_type convention, with ion encoding and adduct metadata where available. Native persistence is the full representation.

`FileHandler` recognizes the native manifest and provides owning load/store overloads. Legacy idXML input can be imported; the adapter defaults to strict export. Compatibility converter APIs warn about information that idXML/mzTab cannot retain. The established PSM Parquet bundle remains a distinct supported input and is not confused with the native format. No reader for discarded experimental versions is provided.

The benchmark in `tools/benchmarks/identification_data` measures owning generation/load, native writing, score-only scans, full scans and streaming filtering in separate Release processes. Run count and match count are configurable. Synthetic scaling and real idXML measurements must be reported separately. Descriptor residency, Arrow buffers and individual large payloads are additional to callback batch targets; byte targets are not a hard process-memory limit.

See the benchmark README for reproducible commands and the accompanying validation report for measured results.

## Checked dataset-wide score schema

A dataset has one ordered PSM score schema and one primary column. Configured runs
must declare exactly the same complete definitions, including direction, scope,
software, parameters and calibration provenance, in the same order. Equal names
alone are insufficient. Protein and inference scores remain independent.

Primary values are required on every match. Supplementary values may be null,
including when an entire run lacks values for a declared supplementary column.
Empty runs without score definitions are construction placeholders and remain
unconfigured through a native roundtrip; their zero-row slices use the shared
physical table schema.

`getScoreDefinitions()` and `getPrimaryScoreDefinition()` check this contract.
`setPrimaryScore(definition)` checks the shared layout and candidate coverage in
every configured run before changing any selection. It can repair differing
primary selections in otherwise compatible runs; failure changes no selections.
`addRun`, `replaceRun`, legacy import, native read/write, and inference reject
incompatible datasets. Native descriptor checks also cover inspection, streaming
scans and filtering before rows are delivered or output is created. Mutable C++
run access supports incremental edits; call `validate()` after those edits.

Normalize incompatible inputs explicitly before combining them. The writer does
not discover a union of scores, reorder columns, demote declared scores to metadata,
or create alternative match files. Additional engine-specific values can remain
in typed match metadata. The legacy importer promotes the primary score and retains
supplementary legacy values in metadata.

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
