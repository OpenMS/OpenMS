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

- A **run** groups identifications with one configuration (search settings in a
  `ProteinIdentification`), one molecule kind (peptide, oligonucleotide or compound) and the
  dataset's score schema. It owns its sources, which record exact file paths.
- A **query** (`Identification`) is an observation, e.g. a spectrum or a feature, with zero or
  more candidate **matches**. Queries without candidates are valid.
- A match owns its molecular representation (a string with an explicit `Encoding`: AASequence or
  NASequence notation, SMILES, InChI or a database identifier), charge, optional adduct, parent
  evidence, peak annotations, metadata and dense score values.
- An **inference result** stores protein and group values once for any number of contributing runs,
  with run-level provenance (`InferenceInput`). It does not store per-PSM assignments.

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
with `identification:catalog=true` in their processing metadata) declare no scores.

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
| `manifest.json` | Format and schema version, run descriptors (UUID, sources, score definitions, metadata descriptors, ID counters) and table slices |
| `queries.parquet` | Queries of all runs |
| `matches.parquet` | Matches with one column per score definition |
| `parents.parquet` | Optional parent catalogs |
| `inputs.parquet`, `proteins.parquet`, `groups.parquet` | Inference provenance, protein hits and groups (members stored inline) |

Runs and inference results share physical tables and row groups. The manifest records each slice
(start row, row count), and the last column of each table records which run or inference result owns
a row: `run_uuid` in `queries`, `matches` and `parents`, and `inference_identifier` (unique within a
dataset) in `inputs`, `proteins` and `groups`. Readers check every row against its slice, so
`SELECT ... FROM 'x.idparquet/matches.parquet' WHERE run_uuid = '...'` needs no manifest lookup.
Score columns follow the dataset schema; runs without scores leave them empty. Matches follow their
query, so readers stream a run without a query-to-match map.

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

The four-table `.idparquet` bundle of OpenMS 3.6 (psms, proteins, protein_groups, search_params)
remains supported as a different layout of the same file type; `FileHandler` recognizes native
bundles by their manifest.

## Feature and consensus maps

`FeatureMap` and `ConsensusMap` own an `IdentificationData` value (`getIdentificationData()`).
Features link to it with `addIDQuery()`, `addIDMatch()` and `setPrimaryID()`, independently of
the legacy `PeptideIdentification`s they may also carry.

`.featureparquet` and `.consensusparquet` bundles store these parts in addition to their tables:

- `identifications/`: the map's identification data as a native bundle;
- `identification_links.parquet`: one row per link (`primary`, `query` or `match`), keyed by
  feature unique ID like `psms.parquet`. Columns: `feature_unique_id`, `link`, `run_uuid`,
  `record_id` (query or match ID), and `encoding`/`representation` for primary molecules.

Every link must resolve and linked features need distinct valid unique IDs; otherwise the export
throws `Exception::InvalidValue`. Map bundles are also written to a temporary sibling directory; an
existing bundle of the same kind or an empty directory is replaced, any other existing path is left
alone. Maps without native identification data write neither part, and OpenMS 3.6 reads the
extended bundles while ignoring the new parts.

AccurateMassSearch writes its ID-format annotations (`-out_annotation *.featureparquet`) this way;
unmatched masses are queries without candidates.

## Conversion to and from the established classes

`IdentificationDataAdapter::fromLegacy` imports peptide/protein identifications; `toLegacy`
exports them. Export is strict by default: information the established classes cannot represent
(e.g. explicit selected candidates, compound or oligonucleotide runs) is rejected.
`LossPolicy::ALLOW` returns a loss report instead.

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
- In pyOpenMS, IDs and references compare and hash by value. Getters return copies; assign edited
  nested values back to the owning record. Runs are edited through `IdentificationData.run_view()`
  (also returned by `addRun`), and the identification data of feature and consensus maps through
  `identification_data_view()`; a copied run is never written back.

## Limits

Owning loads keep every record in memory; datasets larger than memory are processed with `scan` and
`filter` or one run at a time. Inference over many runs may still need substantial memory. A bundle
is rewritten as a whole: there are no partial updates, and changing scores means writing a new
bundle. Arbitrary non-UTF-8 strings are rejected rather than converted.
