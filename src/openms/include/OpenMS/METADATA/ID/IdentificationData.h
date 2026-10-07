// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SearchParameters.h>
#include <array>
#include <atomic>
#include <compare>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
class IdentificationDataFile;

/**
  @brief Owning identification values with one ordered score schema per dataset.

  A run owns sources, observations and candidate matches. Editing or filtering a run
  does not traverse inference provenance. IDs survive copies, filtering and native
  persistence; a new independent run receives a new UUID, and only the run inside the
  dataset allocates new IDs: runs are edited in place through getRun(), and a copied run
  cannot be written back. Views are invalidated by structural edits; match references may
  also be invalidated by match replacement or transformation. Retain IDs across edits. Instances are not safe for concurrent
  mutation. Concurrent const lookups are safe; the first lookup builds a lazy index,
  which Run::prepareLookupIndexes() can build up front.
  @ingroup Metadata
*/
class OPENMS_DLLAPI IdentificationData
{
public:
  enum class MoleculeKind
  {
    PEPTIDE,
    OLIGONUCLEOTIDE,
    COMPOUND
  };
  enum class Encoding
  {
    AA_SEQUENCE,
    NA_SEQUENCE,
    SMILES,
    INCHI,
    DATABASE_ID
  };
  enum class TargetDecoy
  {
    UNKNOWN,
    TARGET,
    DECOY,
    BOTH
  };
  enum class ScoreScope
  {
    MATCH,
    PEPTIDE,
    PROTEIN,
    PROTEIN_GROUP,
    OTHER
  };
  enum class InferencePolicy
  {
    PRESERVE,
    DISCARD
  };

  struct OPENMS_DLLAPI QueryId
  {
    UInt64 value = 0;
    auto operator<=>(const QueryId&) const = default;
  };
  struct OPENMS_DLLAPI MatchId
  {
    UInt64 value = 0;
    auto operator<=>(const MatchId&) const = default;
  };
  struct OPENMS_DLLAPI QueryReference
  {
    std::string run_uuid;
    QueryId query;
    auto operator<=>(const QueryReference&) const = default;
  };
  /// Value identity for a molecule annotation; no reference into a dataset.
  struct OPENMS_DLLAPI MoleculeIdentity
  {
    Encoding encoding = Encoding::AA_SEQUENCE;
    std::string representation;
    auto operator<=>(const MoleculeIdentity&) const = default;
  };
  /// Stable association to an owned match; copying a map preserves this value.
  struct OPENMS_DLLAPI MatchReference
  {
    std::string run_uuid;
    MatchId match;
    auto operator<=>(const MatchReference&) const = default;
  };
  struct OPENMS_DLLAPI ScoreId
  {
    UInt32 value = 0;
    UInt64 owner = 0; ///< Runtime schema guard; not persistent identity.
    auto operator<=>(const ScoreId&) const = default;
  };
  struct OPENMS_DLLAPI SourceId
  {
    UInt32 value = 0;
    UInt64 owner = 0;
    auto operator<=>(const SourceId&) const = default;
  };
  /// Index of a database in the databases of its run (Run::getDatabases()).
  struct OPENMS_DLLAPI DatabaseId
  {
    UInt32 value = 0;
    auto operator<=>(const DatabaseId&) const = default;
  };
  /// An accession together with the name or path of its database; comparable across runs and datasets.
  struct OPENMS_DLLAPI QualifiedAccession
  {
    std::string database;
    std::string accession;
    auto operator<=>(const QualifiedAccession&) const = default;
  };
  struct OPENMS_DLLAPI ScoreDefinition
  {
    std::string name;
    std::string accession;
    bool higher_better = true;
    ScoreScope scope = ScoreScope::MATCH;
    std::string software;
    std::string software_version;
    MetaInfoInterface parameters;
    std::string calibration;
    std::string aggregation;
    bool operator==(const ScoreDefinition&) const = default;
  };
  /**
    @brief How a run was produced: the software, its search settings and further processing metadata

    The files of a run are its sources and its databases are its own records, so the settings name neither:
    the metadata must not list files in 'spectra_data' (an empty list, as legacy runs without files have, may
    stay), and the database fields of @p search (db, db_version, taxonomy) stay empty. The raw files behind the sources ('spectra_data_raw') and processing history
    (e.g. 'alignment:*') are metadata.
  */
  struct OPENMS_DLLAPI RunSettings : MetaInfoInterface
  {
    std::string software;         ///< Search engine or tool that produced the run
    std::string software_version;
    DateTime date;
    SearchParameters search;
    bool operator==(const RunSettings&) const = default;
  };
  /// One file that identifications come from (e.g. an mzML file, or the FASTA of a digest catalog).
  /// An empty @p path stands for a file that is not known.
  struct OPENMS_DLLAPI SourceFile : MetaInfoInterface
  {
    std::string identifier;
    std::string path;
    bool operator==(const SourceFile&) const = default;
  };
  /// A sequence database that a run's matches refer to, e.g. a FASTA file (cf. mzIdentML SearchDatabase).
  struct OPENMS_DLLAPI Database : MetaInfoInterface
  {
    std::string path; ///< Path or name of the database, as the search engine reports it
    std::string version;
    std::string taxonomy;
    bool operator==(const Database&) const = default;
  };
  /// An entry of a database: a protein in peptide runs, a nucleic acid in oligonucleotide runs (cf. mzIdentML DBSequence).
  struct OPENMS_DLLAPI DatabaseSequence : MetaInfoInterface
  {
    DatabaseSequence() = default;
    /// An entry of @p database, e.g. <tt>{db, "P02769|ALBU_BOVIN", TargetDecoy::TARGET}</tt>
    DatabaseSequence(DatabaseId database, std::string accession, TargetDecoy target_decoy = TargetDecoy::UNKNOWN) :
        database(database), accession(std::move(accession)), target_decoy(target_decoy)
    {}

    DatabaseId database;
    std::string accession;
    TargetDecoy target_decoy = TargetDecoy::UNKNOWN;
    std::string sequence;
    std::string description;
    bool operator==(const DatabaseSequence&) const = default;
  };
  /// Where a match occurs in a database sequence of its run (cf. mzIdentML PeptideEvidence); an empty accession is an unknown entry.
  struct OPENMS_DLLAPI SequenceEvidence
  {
    DatabaseId database;
    std::string accession;
    std::optional<UInt64> start;
    std::optional<UInt64> end;
    std::string before;
    std::string after;
    bool operator==(const SequenceEvidence&) const = default;
  };
  struct OPENMS_DLLAPI Observation : MetaInfoInterface
  {
    std::string data_id;
    std::optional<double> rt;
    std::optional<double> mz;
    bool operator==(const Observation&) const = default;
  };
  struct OPENMS_DLLAPI MatchData : MetaInfoInterface
  {
    std::string representation;
    Encoding encoding = Encoding::AA_SEQUENCE;
    Int charge = 0;
    std::optional<double> calculated_mz;
    TargetDecoy target_decoy = TargetDecoy::UNKNOWN;
    std::string name;
    std::optional<std::string> formula;
    std::vector<QualifiedAccession> identifiers;
    std::optional<AdductInfo> adduct;
    std::vector<SequenceEvidence> sequence_evidence;
    std::vector<PeptideHit::PeakAnnotation> peak_annotations;
    bool operator==(const MatchData&) const = default;
  };

  class Run;
  class ScoreView;
  class OPENMS_DLLAPI Match : public MatchData
  {
  public:
    MatchId getId() const
    { return id_; }
    const MatchData& getData() const
    { return *this; }
    std::vector<std::optional<double>> getScores() const;
    /// Dense storage; NaN means a missing score, never a numeric score value.
    const std::vector<double>& getScoreValues() const
    { return scores_; }

  private:
    friend class Run;
    friend class ScoreView;
    MatchId id_;
    UInt64 schema_token_ = 0;
    std::vector<double> scores_;
  };
  class OPENMS_DLLAPI Identification : public Observation
  {
  public:
    QueryId getId() const
    { return id_; }
    const Observation& getObservation() const
    { return *this; }
    const std::vector<Match>& getMatches() const
    { return matches_; }
    std::optional<MatchId> getSelectedMatch() const
    { return selected_; }

  private:
    friend class Run;
    QueryId id_;
    std::vector<Match> matches_;
    std::optional<MatchId> selected_;
  };
  /// One file of a run and the identifications made from it. The sources of a run, in order, are
  /// its file list: a file may appear more than once, and a file without identifications keeps its source.
  struct OPENMS_DLLAPI Source
  {
    SourceId id;
    SourceFile file;
    std::vector<Identification> identifications;
  };
  /// Bound score access avoids repeated definition lookup. Rejects foreign schemas.
  class OPENMS_DLLAPI ScoreView
  {
  public:
    std::optional<double> operator()(const Match& match) const;
    const ScoreDefinition& getDefinition() const
    { return definition_; }

  private:
    friend class Run;
    ScoreId id_;
    UInt64 schema_token_ = 0;
    ScoreDefinition definition_;
  };

  class OPENMS_DLLAPI Run
  {
  public:
    explicit Run(std::string identifier = {}, MoleculeKind kind = MoleculeKind::PEPTIDE);
    Run(const Run&);
    Run(Run&&);
    /// Not assignable: a run inside a dataset is edited in place (getRun), so a stale copy can never
    /// replace it and reuse IDs that the live run has allocated meanwhile.
    Run& operator=(const Run&) = delete;
    Run& operator=(Run&&) = delete;
    ~Run() = default;
    const std::string& getIdentifier() const
    { return identifier_; }
    const std::string& getUuid() const
    { return uuid_; }
    MoleculeKind getMoleculeKind() const
    { return kind_; }
    const RunSettings& getSettings() const
    { return *settings_; }
    /// @throw Exception::InvalidValue if the metadata of @p settings lists files in 'spectra_data' (the files of a run are its sources)
    void setSettings(const RunSettings& settings);
    /// The databases that the database sequences and the sequence evidence of the run refer to.
    const std::vector<Database>& getDatabases() const
    { return databases_; }
    /// Add a database, or return the ID of an equal one. Databases are never removed.
    DatabaseId addDatabase(const Database& database);
    DatabaseId getDatabaseId(UInt32 index) const;
    const Database& getDatabase(DatabaseId database) const;
    /// The identity of an entry of a database of the run, comparable across runs (database by path).
    QualifiedAccession qualify(DatabaseId database, const std::string& accession) const;
    /// The database entries that the run's matches refer to, if the run has a catalogue of them.
    const std::optional<std::vector<DatabaseSequence>>& getDatabaseSequences() const
    { return sequences_; }
    /// @throw Exception::InvalidValue for an unknown database, an empty accession or a duplicate (database, accession)
    void setDatabaseSequences(std::optional<std::vector<DatabaseSequence>> sequences);
    const std::vector<Source>& getSources() const
    { return sources_; }
    const std::vector<ScoreDefinition>& getScoreDefinitions() const
    { return scores_; }
    /// Append a file to the run's file list. Sources are never removed, so their order stays stable.
    SourceId addSource(const SourceFile& source);
    SourceId getSourceId(UInt32 index) const;
    ScoreId addScore(const ScoreDefinition& definition);
    ScoreId getScoreId(UInt32 index) const;
    ScoreId findScore(const ScoreDefinition& definition) const;
    const ScoreDefinition& getScoreDefinition(ScoreId score) const;
    ScoreView bindScore(ScoreId score) const;
    std::optional<ScoreId> getPrimaryScore() const
    { return primary_; }
    void setPrimaryScore(std::optional<ScoreId> score);
    QueryId addIdentification(SourceId source, const Observation& observation);
    MatchId addMatch(QueryId query, const MatchData& data, const std::vector<std::optional<double>>& scores = {});
    const Identification* findIdentification(QueryId id) const;
    const Match* findMatch(MatchId id) const;
    const Identification& getIdentification(QueryId id) const;
    const Match& getMatch(MatchId id) const;
    /// Resolve the owning observation without a full dataset scan.
    const Identification& getIdentificationForMatch(MatchId id) const;
    std::optional<double> getScore(MatchId match, ScoreId score) const;
    void setScore(MatchId match, ScoreId score, std::optional<double> value);
    void setSelectedMatch(QueryId query, std::optional<MatchId> selected);
    void replaceObservation(QueryId query, const Observation& observation);
    /// Payload-only edits preserve scores and may not change molecular/ion identity.
    void replaceMatch(MatchId match, const MatchData& data);
    /// Replace the hypothesis and explicitly supply all scores for its new identity.
    void replaceMatch(MatchId match, const MatchData& data, const std::vector<std::optional<double>>& scores);
    /// Evaluate all predicates before committing; throwing callbacks leave values unchanged.
    Size filterMatches(const std::function<bool(const Match&)>& keep, bool keep_empty_queries = false);
    Size eraseMatches(const std::function<bool(const Match&)>& remove, bool keep_empty_queries = false);
    /// Remove identifications (with their matches), also ones without matches; returns the number removed.
    /// Evaluates all predicates before committing, like filterMatches().
    Size eraseIdentifications(const std::function<bool(const Identification&)>& remove);
    Size retainBest(ScoreId score, bool keep_ties = true, bool keep_empty_queries = false);
    void transformMatches(const std::function<void(MatchData&)>& transform);
    Size getNumberOfIdentifications() const;
    Size getNumberOfMatches() const;
    UInt64 getNextQueryId() const
    { return next_query_id_; }
    UInt64 getNextMatchId() const
    { return next_match_id_; }
    /// Build the lazy ID lookup indexes up front (const lookups otherwise build them on first use); mutations require exclusive access.
    void prepareLookupIndexes();
    /// Release the spare capacity of the run's containers (sources, queries, candidates), e.g. after importing records one
    /// by one. Like any structural edit, this invalidates record views; IDs and lookup indexes stay valid.
    void shrinkToFit();
    /// Import explicit IDs during construction. After restoration/filtering, historical IDs cannot be reused.
    /// Payloads are owned by value so importers can transfer decoded records without copying.
    QueryId importIdentification(SourceId source, QueryId id, Observation observation);
    MatchId importMatch(QueryId query, MatchId id, MatchData data, const std::vector<std::optional<double>>& scores = {});
    /// Restore persisted UUID and counters; counters must exceed every live ID.
    void restoreIdentity(const std::string& uuid, UInt64 next_query, UInt64 next_match);
    /// Reserve IDs appearing only in retained inference provenance.
    void reserveMatchId(MatchId id);
    void validate() const;
    bool operator==(const Run& other) const;

  private:
    friend class IdentificationData;
    std::string identifier_;
    std::string uuid_;
    MoleculeKind kind_;
    // Standard-library trees inside the settings can allocate when moved on MSVC.
    // Indirection keeps the run's transactional commit nonthrowing.
    std::unique_ptr<RunSettings> settings_ = std::make_unique<RunSettings>();
    std::vector<Database> databases_;
    std::optional<std::vector<DatabaseSequence>> sequences_;
    std::vector<Source> sources_;
    std::vector<ScoreDefinition> scores_;
    std::vector<UInt64> score_owners_;
    std::optional<ScoreId> primary_;
    UInt64 schema_token_;
    UInt64 next_query_id_ = 1;
    UInt64 next_match_id_ = 1;
    mutable bool callback_active_ = false;
    bool import_finalized_ = false;
    Size query_count_ = 0;
    Size match_count_ = 0;
    // Lazy lookup indexes may be built from concurrent const lookups: the flags are published
    // with release/acquire and the build is serialized; mutations still require exclusive access.
    mutable std::atomic<bool> query_index_built_ {false};
    mutable std::atomic<bool> match_index_built_ {false};
    mutable std::mutex index_mutex_;
    mutable std::unordered_map<UInt64, std::array<Size, 2>> query_index_;
    mutable std::unordered_map<UInt64, std::array<Size, 3>> match_index_;
    std::optional<std::array<Size, 2>> last_query_;
    std::optional<std::array<Size, 3>> last_match_;
    void invalidateIndexes_();
    void ensureQueryIndex_() const;
    void ensureMatchIndex_() const;
    void checkMutation_() const;
    void checkScore_(ScoreId score) const;
    void validateMatch_(const MatchData& data, const std::vector<std::optional<double>>& scores) const;
    void validateMatchData_(const MatchData& data) const;
    Identification& query_(QueryId id);
    Match& match_(MatchId id);
    void swapData_(Run& other) noexcept;
  };

  /// An identification of a run with some of its matches, e.g. those that a feature links. Points into the
  /// dataset, so it is invalidated like a record view by structural edits.
  struct OPENMS_DLLAPI QueryMatches
  {
    const Run* run = nullptr;
    const Identification* query = nullptr;
    /// In the order of the identification's matches
    std::vector<const Match*> matches;
    /// The match with the best primary score (the first of equal ones; matches without a value are skipped), or nullptr
    const Match* getBestMatch() const;
  };

  /// Run-level provenance for an inference calculation; no per-match input list is retained.
  struct OPENMS_DLLAPI InferenceInput
  {
    std::string run_identifier;
    std::string run_uuid;
    std::optional<ScoreDefinition> score;
    /// Description of the selection used at calculation time, not an executable filter.
    std::string selection;
    bool operator==(const InferenceInput&) const = default;
  };
  struct OPENMS_DLLAPI InferenceResult
  {
    std::string identifier;
    ProteinIdentification proteins;
    std::optional<ScoreDefinition> protein_score;
    std::optional<ScoreDefinition> group_score;
    std::map<std::string, QualifiedAccession> qualified_accessions;
    std::vector<InferenceInput> inputs;
    bool operator==(const InferenceResult&) const = default;
  };

  IdentificationData() = default;
  IdentificationData(const IdentificationData&);
  IdentificationData(IdentificationData&&);
  IdentificationData& operator=(const IdentificationData&);
  IdentificationData& operator=(IdentificationData&&);
  ~IdentificationData() = default;
  /// Adding runs preserves existing run references. Record views are invalidated by structural edits.
  Run& addRun(const std::string& identifier, MoleculeKind kind = MoleculeKind::PEPTIDE);
  /// Add an owning copy, rejecting duplicate UUIDs or display identifiers.
  Run& addRun(Run run);
  Run& getRun(const std::string& identifier);
  const Run& getRun(const std::string& identifier) const;
  Run* findRunByUuid(const std::string& uuid);
  const Run* findRunByUuid(const std::string& uuid) const;
  const std::deque<Run>& getRuns() const
  { return runs_; }
  const std::vector<InferenceResult>& getInferenceResults() const
  { return inference_; }
  void addInferenceResult(InferenceResult result);
  void clearInferenceResults();
  bool empty() const
  { return runs_.empty() && inference_.empty(); }
  void clear();
  /// Append independent runs atomically; identical shared UUIDs are retained once.
  /// Conflicting values for an existing UUID are rejected; repeated display names receive a numeric suffix.
  /// Existing runs are neither copied nor moved, so references to them stay valid; a rejected merge changes nothing.
  void merge(const IdentificationData& other);
  /**
    @brief Resolve links to identifications and matches (e.g. of a feature)

    Every linked identification appears once, with its linked matches (an identification linked only by itself has none),
    ordered by identification ID, then run: the order of peptide identifications exported from imported ones.

    @throw Exception::MissingInformation if a link refers to an identification or match that does not exist
  */
  std::vector<QueryMatches> resolveLinks(const std::set<QueryReference>& queries, const std::set<MatchReference>& matches) const;
  /**
    @brief The identifications and matches that none of the given links refers to (e.g. those not assigned to features)

    An identification appears with its matches that are not linked, if it has any; an identification without matches
    appears if it is not linked itself. The order is the order of identification IDs (by run, then ID, if runs share IDs),
    the order of exported unassigned peptide identifications.
  */
  std::vector<QueryMatches> getUnlinked(const std::set<QueryReference>& queries, const std::set<MatchReference>& matches) const;
  bool operator==(const IdentificationData& other) const;
  Size filterMatches(const std::function<bool(const Match&)>& keep, InferencePolicy policy, bool keep_empty_queries = false);
  /** Dataset-wide ordered PSM score contract.
      Configured runs have identical complete ScoreDefinitions in identical column
      order and select the same primary column. Primary values are required;
      supplementary values may be missing. Empty runs with no score definitions
      are construction placeholders and are ignored. Run-local ScoreId handles
      remain distinct even though their column indices agree.
      Throws on disagreement. Mutable run edits must be followed by validate();
      import, merge, export and inference boundaries enforce this contract.
  */
  const std::vector<ScoreDefinition>& getScoreDefinitions() const;
  /// Return the common primary definition after checking the complete score contract.
  std::optional<ScoreDefinition> getPrimaryScoreDefinition() const;
  /// Select an existing score in every participating run, checking coverage first.
  /// Failure leaves all primary selections unchanged. Empty unconfigured runs are ignored.
  void setPrimaryScore(const ScoreDefinition& definition);
  void validate() const;
  void swap(IdentificationData& other);

private:
  std::deque<Run> runs_;
  std::vector<InferenceResult> inference_;
  bool callback_active_ = false;
  void checkMutation_() const;
};
} // namespace OpenMS
