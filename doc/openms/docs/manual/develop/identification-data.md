Identification data model and native persistence
================================================

Since OpenMS 3.7, `IdentificationData` is an owning model of search results. It replaces the
reference graph of earlier versions (`IDDataContainer`, `ObservationMatch`, `IdentifiedMolecule`,
`ScoreType` and related classes) and the OpenMS SQLite format (`.oms`). This page describes the
model, its score contract and its Parquet persistence for developers. The established
`PeptideIdentification`/`ProteinIdentification` classes remain the model of most tools;
`IdentificationDataAdapter` converts between both.

## Model

An `IdentificationData` dataset contains analysis runs and independent inference results.

- A **run** groups identifications with its settings (`RunSettings`: software and version, date,
  `SearchParameters` and metadata), one molecule kind (peptide, oligonucleotide or compound) and the
  dataset's score schema.
- The **sources** of a run are its files, in order. A `Source` is one file (`SourceFile`) with the
  queries from it; a file may appear twice, a file without queries keeps its source, and an empty
  path stands for a file that is not known. The position of a source is the index that legacy PSMs
  store as `id_merge_index`, so the native model needs neither that index nor the `spectra_data`
  list; the run settings must not list `spectra_data` (`spectra_data_raw` is settings metadata).
- The **databases** of a run (`Database`: path, version, taxonomy; cf. mzIdentML `SearchDatabase`)
  are its own records, so the search parameters in the settings leave `db`, `db_version` and
  `taxonomy` empty. Sequence evidence and the optional **database sequences** of a run
  (`DatabaseSequence`: accession, target/decoy, sequence, description; cf. mzIdentML `DBSequence`)
  refer to a database by `DatabaseId`, its index in the run. `Run::qualify` turns a database and an
  accession into a `QualifiedAccession` (database path and accession), which compares across runs
  and datasets; inference results use it.
- A **query** (`Identification`) is an observation, e.g. a spectrum or a feature, with zero or
  more candidate **matches**. Queries without candidates are valid.
- A match owns its molecular representation (a string with an explicit `Encoding`: AASequence or
  NASequence notation, SMILES, InChI or a database identifier), charge, sequence evidence
  (`SequenceEvidence`: database, accession, 0-based positions and one-character flanking residues;
  cf. mzIdentML `PeptideEvidence`), peak annotations and metadata. Name, formula, database
  identifiers and adduct (`MoleculeDetails`, mostly of compounds) are stored apart from the match
  (`MatchData::details`, empty for most peptides).
- The scores of a run are stored in its own columns, one per score definition, with a row per match
  (`Run::getScores()`, `Run::getScore()`, or a `ScoreView` from `Run::bindScore()`). A view shares
  ownership of the columns and belongs to one state of them: adding a score, filtering or copying the
  run starts a new state, and the view (and a copy of a match taken before) is then rejected rather
  than reading another match's value. `Run::getRevision()` counts the edits of a run in a process.
- An **inference result** stores protein and group values once for any number of contributing runs,
  with run-level provenance (`InferenceInput`) and its protein and group score definitions. It does
  not store per-PSM assignments.

Runs have a UUID; queries and matches have stable `UInt64` IDs within their run. IDs survive
copies, filtering and persistence and are never reused. Features refer to identifications by value
(`QueryReference`, `MatchReference`: run UUID plus record ID), so copying a map needs no reference
translation. `IdentificationData::merge` appends runs; references to existing runs stay valid, and a
rejected merge changes nothing. Only the run inside a dataset allocates new IDs: runs are edited in
place through `getRun()`, and `Run` is not assignable, so a stale copy cannot replace a run and
reuse IDs that the run allocated meanwhile.

Concurrent const lookups by ID are safe; the first lookup builds an index per run
(`Run::prepareLookupIndexes()` builds it up front). Mutations need exclusive access.

## Score contract

A dataset has one ordered PSM score schema and one primary score column. Every run with matches
declares exactly the same `ScoreDefinition`s in the same order: name, direction, scope, producing
software and version, parameters and calibration. Equal names alone are not enough. Primary values
are required on every match; supplementary values may be missing. Scoreless sequence catalogs (runs
with `identification:catalog=true` in their settings metadata) declare no scores.

Data that does not satisfy the contract is rejected, and the error names the conflicting
definitions. The model does not split runs or demote scores on its own. To combine searches from
different engines or engine versions, normalize the scores first, for example:

- rescore each search with PercolatorAdapter (q-values/PEPs produced by Percolator), or
- compute posterior error probabilities with IDPosteriorErrorProbability, or
- select a common score with IDScoreSwitcher.

A derived score belongs to the tool that computed it, not to the search engine.
IDPosteriorErrorProbability and FalseDiscoveryRate record themselves with
`ProteinIdentification::setScoreSoftware()` (search parameters `ScoreSoftware:<score type>`),
because downstream tools still need the original search engine. When importing idXML,
`IdentificationDataAdapter` takes the score's producer from that record and falls back to the
search engine. PercolatorAdapter and ConsensusID record themselves as search engine and clear such
records. As a result, PEPs computed by the same IDPosteriorErrorProbability version for Comet and
MS-GF+ searches share one definition, while raw engine scores do not.

## Native identification bundles

`IdentificationDataFile` stores a dataset as a directory (`.idparquet`) with a JSON manifest and
shared Parquet tables:

| Path | Contents |
| --- | --- |
| `manifest.json` | Format and schema version, score column names, run descriptors (UUID, settings, sources, databases, score definitions, metadata descriptors, ID counters) and table slices |
| `queries.parquet` | Queries of all runs |
| `matches.parquet` | Matches with one column per score definition (`score_pep`, `score_q_value`, …) |
| `database_sequences.parquet` | Optional database sequences (`database` is the index of a database of the run) |
| `inputs.parquet`, `proteins.parquet`, `groups.parquet` | Inference provenance, protein hits and groups (members stored inline) |

Runs and inference results share physical tables and row groups. The manifest records each slice
(start row, row count), and the last column of each table records which run or inference result owns
a row: `run_uuid` in `queries`, `matches` and `database_sequences`, and `inference_identifier` (unique within a
dataset) in `inputs`, `proteins` and `groups`. Readers check every row against its slice, so
`SELECT ... FROM 'x.idparquet/matches.parquet' WHERE run_uuid = '...'` needs no manifest lookup.
Score columns follow the dataset schema; runs without scores leave them empty. Matches follow their
query, so readers stream a run without a query-to-match map.

Score columns are named after their definitions: `score_` plus the name in lower case, with every
character other than ASCII letters and digits replaced by `_` (`q-value` → `score_q_value`,
`MS:1002252` → `score_ms_1002252`). Definitions that share a name get their producer appended
(`score_pep_idposteriorerrorprobability`, `score_pep_percolator`), and a number if that is not
enough. Such names need no quoting in SQL. DuckDB in particular reads an unquoted `MS:1002252` as an
alias followed by a number, and treats names that differ only in case as the same. The manifest
records the names (`score_columns`, in score schema order), readers locate scores by them, and
`IdentificationDataFile::inspect` returns them per run. Because equal definitions get equal names,
DuckDB's `union_by_name` combines the same score across bundles.

Run UUIDs are random version-4 UUIDs written as 36-character lowercase strings
(`4ee928b7-e4ed-4925-8899-133a43050310`). Every table that carries them, including
`identification_links.parquet` of map bundles, stores them as dictionary-encoded strings
(`dictionary<int32, string>`, one entry per run) with the Arrow schema stored in the file, so readers
get the full string at about 4 bytes per row; plain string columns are accepted on input. Inference
identifiers are stored the same way.

Writes go to a temporary sibling directory that is renamed into place when every table is closed.
An existing destination is rejected unless `Options::replace_existing` is set, which replaces an
existing native bundle (and nothing else). `FileHandler::storeIdentifications(filename,
IdentificationData)` sets it, so a tool rerun overwrites its previous output. Owning loads validate
into a temporary dataset and leave the destination unchanged on failure.

Besides owning load and store, the API offers projected streaming scans (`scan`), single-run loading
(`loadRun`), descriptor inspection (`inspect`) and streaming subset export into a new bundle
(`filter`), which keeps IDs and either preserves or discards inference results. `Options::threads`
enables a private Arrow thread pool per operation; the default is serial.

`.idparquet` is the native format only. The `FileHandler` overloads for the established classes
import and export through `IdentificationDataAdapter` (below), so every tool that reads or writes
`.idparquet` uses native bundles. The four-table bundle of OpenMS 3.6 (psms, proteins,
protein_groups, search_params) is not read: there is no importer, and loading one reports that it
must be converted to idXML with OpenMS 3.6. OpenMS 3.6 does not read native bundles either.

## Feature and consensus maps

`FeatureMap` and `ConsensusMap` own an `IdentificationData` value (`getIdentificationData()`).
Features link to it with `addIDQuery()`, `addIDMatch()` and `setPrimaryID()`, independently of
the legacy `PeptideIdentification`s they may also carry.

`.featureparquet` and `.consensusparquet` bundles store these parts in addition to their tables:

- `identifications/`: the map's identification data as a native bundle;
- `identification_links.parquet`: one row per link (`primary`, `query` or `match`), keyed by
  feature unique ID like `psms.parquet`. Columns: `feature_unique_id`, `link`, `run_uuid`,
  `record_id` (query or match ID), and `encoding`/`representation` for primary molecules.
  `feature_unique_id` is an `int64` holding the bits of the unsigned unique ID, like `unique_id` in
  the feature tables and `feature_unique_id` in `psms.parquet`, so IDs of 2^63 and above appear as
  negative numbers and the tables join without a cast.

Links join the feature table and the identification tables without the manifest. For example, the
best accurate-mass annotations per feature with DuckDB:

```sql
SELECT f.unique_id, f.rt, f.mz, m.name, m.adduct.name AS adduct, m.score_masserrorppmscore AS ppm
FROM 'x.featureparquet/features.parquet' f
JOIN 'x.featureparquet/identification_links.parquet' l
  ON l.feature_unique_id = f.unique_id AND l.link = 'match'
JOIN 'x.featureparquet/identifications/matches.parquet' m
  ON m.run_uuid = l.run_uuid AND m.match_id = l.record_id
WHERE abs(m.score_masserrorppmscore) < 2
ORDER BY f.unique_id, abs(m.score_masserrorppmscore);
```

Query links join `identifications/queries.parquet` on `run_uuid` and `query_id` the same way.

Every link must resolve and linked features need distinct valid unique IDs; otherwise the export
throws `Exception::InvalidValue`. Map bundles are also written to a temporary sibling directory; an
existing bundle of the same kind or an empty directory is replaced, any other existing path is left
alone. Maps without native identification data write neither part.

AccurateMassSearch writes its ID-format annotations (`-out_annotation *.featureparquet`) this way;
unmatched masses are queries without candidates.

featureXML and consensusXML hold legacy identifications only. A map with native identification data
and no legacy identifications is written with the native ones converted
(`IdentificationDataConverter::exportFeatureIDs`/`exportConsensusIDs`) instead of losing them; a map
with legacy identifications is written as before, without its native data.

## Conversion to and from the established classes

`IdentificationDataAdapter::fromLegacy` imports peptide/protein identifications; `toLegacy`
exports them. Export is strict by default: information the established classes cannot represent
(e.g. explicit selected candidates, compound or oligonucleotide runs, a run with more than one
database) is rejected. `LossPolicy::ALLOW` returns a loss report instead.

A legacy protein run maps to a run as follows: search engine, version, date, search parameters and
meta values become the settings (`settingsFromLegacy`); `db`, `db_version` and `taxonomy` become a
database of the run (`databaseFromLegacy`); `spectra_data` becomes the sources, and each PSM goes to
the source its `id_merge_index` names, else to the only file, else to a source without a path
(`addLegacySources`, `legacySource`); the protein hits become database sequences. The protein run is
kept as the inference result `legacy:<run>` only if export cannot rebuild it exactly from the run:
search engine output, which lists the proteins of its matches without scores, is not kept, while
protein scores, groups or coverage are. A protein score type other than the PSM score (often empty,
or the search engine score after rescoring) is then kept in the settings metadata
(`identification:legacy_protein_score_type`, `identification:legacy_protein_higher_score_better`). Export
reverses this (`settingsToLegacy`, `legacyFiles`) and writes `id_merge_index` only
for runs with several files.

A round trip keeps the order of the legacy values. Import creates the runs in the order of the
protein runs, and numbers the queries of all runs in the order of the peptide identifications.
Export writes the protein runs in run order and the peptide identifications by query ID across all
runs, which restores the legacy order; if two runs share a query ID (runs created independently, e.g.
merged datasets), it writes them run by run instead, each by query ID. The golden outputs of
IDMerger, PercolatorAdapter, ProteomicsLFQ and IsobaricWorkflow survive legacy → native → legacy
unchanged (`IdentificationDataRoundTrip_test`). Legacy meta values that a field holds are not stored twice: the spectrum
reference is the query's `data_id`, `target_decoy` of peptide and protein hits is the
`target_decoy` field of matches and database sequences, and the description of a protein hit is that
of its database sequence. Export writes them back; a meta value that
export would write differently (another spelling or value type) is kept as metadata.

An inference result becomes one protein run. Inference over several runs is exported as the
legacy model represents it: one merged protein run (as IDMerger creates it) that lists the files of
all input runs, with `id_merge_index` on each PSM. A run joins only if its search engine and
settings are mergeable by the legacy rules (`SearchParameters::mergeable`); settings that are
mergeable but not identical are reported as a loss, and a run that cannot join is exported as its
own protein run without the inference result, never under another run's settings. The score
definitions of the inference result (protein, group and input score) are stored as
`identification:inference:*` metadata of the protein run, and import restores them.

`IdentificationDataConverter` keeps the RNA and compound conventions of idXML and mzTab, and
converts feature annotations (`importFeatureIDs`, `exportFeatureIDs`).

`FileHandler::loadIdentifications(filename, IdentificationData&)` reads native bundles directly and
imports other identification formats through the adapter.
`FileHandler::storeIdentifications(filename, const IdentificationData&)` writes `.idparquet` natively
and other formats through strict export.

TOPP tools using the owning model: IDMerger and MapAlignerIdentification merge or align native
bundles when all inputs are native; NucleicAcidSearchEngine writes its digest and results as native
bundles (`-digest_out`, `-db_out`) and reads digests (`-digest`); AccurateMassSearch as above.

## Migration

- `.oms` files are no longer supported. Convert them with OpenMS 3.6: IDFileConverter for
  identifications (e.g. to idXML) and FileConverter for feature maps (to featureXML).
- Code using the reference graph needs to move to IDs: `ObservationMatchRef` becomes
  `MatchReference`, iteration over `getObservationMatches()` becomes iteration over runs, sources,
  queries and matches, and score lookups go through `ScoreId` or a bound `ScoreView`.
  `BaseFeature::updateIDReferences()` is no longer needed.
- In pyOpenMS, all value types compare by value (`==`), take their fields as constructor keywords and
  print them. IDs (including `DatabaseId`), references, `MoleculeIdentity` and `QualifiedAccession`
  also hash by value; mutable records are unhashable. Getters return copies; assign edited nested
  values back to the owning record. Runs are edited through `IdentificationData.run_view()`
  (also returned by `addRun`), and the identification data of feature and consensus maps through
  `identification_data_view()`; a copied run is never written back.

## Limits

A match takes 136 bytes on 64-bit Linux (a `PeptideHit` 120) and a sequence evidence 64 (a
`PeptideEvidence` 48); scores take 8 bytes per value in the columns of their run. Imported protein
lists are stored once: the database sequences hold sequence, description and metadata, and the
`legacy:<run>` inference result only the inference values of the proteins.

Owning loads keep every record in memory; datasets larger than memory are processed with `scan` and
`filter` or one run at a time. Loading a bundle or importing legacy identifications ends with
`Run::shrinkToFit()`, which releases the spare capacity of containers filled one record at a time;
code that builds large runs record by record can call it as well. Inference over many runs may still need substantial memory. A bundle
is rewritten as a whole: there are no partial updates, and changing scores means writing a new
bundle. Arbitrary non-UTF-8 strings are rejected rather than converted.
