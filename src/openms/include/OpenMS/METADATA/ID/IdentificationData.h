// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <array>
#include <compare>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
class IdentificationDataFile;

/**
  @brief Owning identification values with a homogeneous score contract per analysis run.

  A run owns sources, observations and candidate matches. Editing or filtering a run
  does not traverse inference provenance. IDs survive copies, filtering and native
  persistence; a new independent run receives a new UUID. Views are invalidated by
  structural edits. Instances are not safe for concurrent mutation. Prepare lazy lookup
    indexes with Run::prepareLookupIndexes() before parallel const ID lookups.
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
  struct OPENMS_DLLAPI SourceFile : MetaInfoInterface
  {
    std::string identifier;
    std::string path;
    std::vector<std::string> primary_files;
    bool operator==(const SourceFile&) const = default;
  };
  struct OPENMS_DLLAPI ParentEvidence
  {
    QualifiedAccession parent;
    std::optional<UInt64> start;
    std::optional<UInt64> end;
    std::string before;
    std::string after;
    bool operator==(const ParentEvidence&) const = default;
  };
  struct OPENMS_DLLAPI ParentRecord : MetaInfoInterface
  {
    QualifiedAccession identity;
    TargetDecoy target_decoy = TargetDecoy::UNKNOWN;
    std::string sequence;
    std::string description;
  };
  struct OPENMS_DLLAPI Observation : MetaInfoInterface
  {
    std::string data_id;
    std::optional<double> rt;
    std::optional<double> mz;
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
    std::vector<ParentEvidence> parent_evidence;
    std::vector<PeptideHit::PeakAnnotation> peak_annotations;
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
  struct OPENMS_DLLAPI SourceBlock
  {
    SourceId id;
    SourceFile source;
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
    Run& operator=(const Run&);
    Run& operator=(Run&&);
    ~Run() = default;
    const std::string& getIdentifier() const
    { return identifier_; }
    const std::string& getUuid() const
    { return uuid_; }
    MoleculeKind getMoleculeKind() const
    { return kind_; }
    const ProteinIdentification& getProcessingMetadata() const
    { return processing_; }
    void setProcessingMetadata(const ProteinIdentification& metadata);
    const std::optional<std::vector<ParentRecord>>& getParents() const
    { return parents_; }
    void setParents(std::optional<std::vector<ParentRecord>> parents);
    const std::vector<SourceBlock>& getSourceBlocks() const
    { return sources_; }
    const std::vector<ScoreDefinition>& getScoreDefinitions() const
    { return scores_; }
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
    Size retainBest(ScoreId score, bool keep_ties = true, bool keep_empty_queries = false);
    void transformMatches(const std::function<void(MatchData&)>& transform);
    Size getNumberOfIdentifications() const;
    Size getNumberOfMatches() const;
    UInt64 getNextQueryId() const
    { return next_query_id_; }
    UInt64 getNextMatchId() const
    { return next_match_id_; }
    /// Prepare optional indexes before parallel const lookups; mutations require exclusive access.
    void prepareLookupIndexes();
    /// Import explicit IDs during construction. After restoration/filtering, historical IDs cannot be reused.
    QueryId importIdentification(SourceId source, QueryId id, const Observation& observation);
    MatchId importMatch(QueryId query, MatchId id, const MatchData& data, const std::vector<std::optional<double>>& scores = {});
    /// Restore persisted UUID and counters; counters must exceed every live ID.
    void restoreIdentity(const std::string& uuid, UInt64 next_query, UInt64 next_match);
    /// Reserve IDs appearing only in retained inference provenance.
    void reserveMatchId(MatchId id);
    void validate() const;

  private:
    friend class IdentificationData;
    std::string identifier_;
    std::string uuid_;
    MoleculeKind kind_;
    ProteinIdentification processing_;
    std::optional<std::vector<ParentRecord>> parents_;
    std::vector<SourceBlock> sources_;
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
    mutable bool query_index_built_ = false;
    mutable bool match_index_built_ = false;
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
    Identification& query_(QueryId id);
    Match& match_(MatchId id);
    void swapData_(Run& other) noexcept;
  };

  /// Run-level provenance for an inference calculation; no per-match input list is retained.
  struct OPENMS_DLLAPI InferenceInput
  {
    std::string run_identifier;
    std::string run_uuid;
    std::optional<ScoreDefinition> score;
    /// Description of the selection used at calculation time, not an executable filter.
    std::string selection;
  };
  struct OPENMS_DLLAPI InferenceResult
  {
    std::string identifier;
    ProteinIdentification proteins;
    std::optional<ScoreDefinition> parent_score;
    std::optional<ScoreDefinition> group_score;
    std::map<std::string, QualifiedAccession> qualified_accessions;
    std::vector<InferenceInput> inputs;
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
  /// Replace an existing run of identical persistent identity, preserving inference.
  void replaceRun(const Run& run);
  Run& getRun(const std::string& identifier);
  const Run& getRun(const std::string& identifier) const;
  Run* findRunByUuid(const std::string& uuid);
  const Run* findRunByUuid(const std::string& uuid) const;
  const std::deque<Run>& getRuns() const
  { return runs_; }
  const std::vector<InferenceResult>& getInferenceResults() const
  { return inference_; }
  void addInferenceResult(const InferenceResult& result);
  void clearInferenceResults();
  Size filterMatches(const std::function<bool(const Match&)>& keep, InferencePolicy policy, bool keep_empty_queries = false);
  /** Common primary PSM score contract.
      All runs with matches or a selected primary score must agree on the complete
      ScoreDefinition (including orientation and provenance), or all be unscored.
      Empty runs without a primary score are construction placeholders and are ignored.
      Supplementary score definitions and local score IDs may differ between runs.
      Throws on disagreement. Mutable run edits must be followed by validate();
      import, replacement, export and inference boundaries enforce this contract.
  */
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
