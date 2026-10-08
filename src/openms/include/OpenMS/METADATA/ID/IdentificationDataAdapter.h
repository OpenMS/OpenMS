// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <functional>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

namespace OpenMS
{
class FeatureMap;
class ConsensusMap;

/**
  @brief Explicit conversions between owning identifications and established peptide/protein values.

  Imports reject heterogeneous primary PSM score definitions; normalize them first. Export
  rejects unrepresentable information by default. ALLOW returns a loss report instead of
  silently discarding information. No inference freshness is inferred from current matches.
  @ingroup Metadata
*/
class OPENMS_DLLAPI IdentificationDataAdapter
{
public:
  enum class LossPolicy
  {
    STRICT,
    ALLOW
  };
  enum class MissingLinkPolicy
  {
    REJECT,
    PRUNE
  };
  struct OPENMS_DLLAPI ExportOptions
  {
    LossPolicy loss_policy = LossPolicy::STRICT;
    bool include_inference = true;
    std::optional<std::string> inference_result;
  };
  using QueryReference = IdentificationData::QueryReference;
  struct OPENMS_DLLAPI ImportResult
  {
    IdentificationData data;
    std::vector<QueryReference> queries; ///< In the original peptide-identification order.
  };
  struct OPENMS_DLLAPI LegacyResult
  {
    std::vector<ProteinIdentification> proteins;
    PeptideIdentificationList peptides;
    std::vector<QueryReference> queries; ///< Parallel to peptides.
    std::vector<std::string> losses;
  };
  struct OPENMS_DLLAPI FeatureAssociation
  {
    bool unassigned = false;
    std::vector<Size> feature_path; ///< Top-level index followed by subordinate indices.
    QueryReference query;
    std::vector<IdentificationData::MatchId> matches;
  };
  struct OPENMS_DLLAPI FeatureImportResult
  {
    IdentificationData data;
    std::vector<FeatureAssociation> associations;
  };

  /// Import with exact observation mappings, including identifications with no hits.
  static ImportResult importLegacy(const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides);
  /// Convenience import without the optional legacy-order mapping.
  static IdentificationData fromLegacy(const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides);
  /// Strict export with automatic selection only when an inference result is unambiguous.
  static LegacyResult toLegacy(const IdentificationData& data);
  static LegacyResult toLegacy(const IdentificationData& data, const ExportOptions& options);
  /// Explicit chemical materialization using the run's retained custom modification definitions.
  static PeptideHit materializePeptide(const IdentificationData::Run& run, const IdentificationData::Match& match, IdentificationData::ScoreId score);
  /// Collect live associations while retaining map values outside the identification model.
  static FeatureImportResult fromFeatureMap(const FeatureMap& map);
  static FeatureImportResult fromConsensusMap(const ConsensusMap& map);
  /// Prune deleted candidate links, or throw without changing any associations.
  static Size reconcileAssociations(const IdentificationData& data, std::vector<FeatureAssociation>& associations, MissingLinkPolicy policy);
  /// Replace identification annotations only; retain measured features, channels and intensities.
  static std::vector<std::string> applyToFeatureMap(const IdentificationData& data,
                                                    const std::vector<FeatureAssociation>& associations,
                                                    FeatureMap& map,
                                                    const ExportOptions& options,
                                                    MissingLinkPolicy policy);
  static std::vector<std::string> applyToConsensusMap(const IdentificationData& data,
                                                      const std::vector<FeatureAssociation>& associations,
                                                      ConsensusMap& map,
                                                      const ExportOptions& options,
                                                      MissingLinkPolicy policy);

  /**
    @brief Run settings from a legacy protein run: search engine and version, date, search parameters and metadata

    The protein values of the legacy run belong to an inference result, its identifier is the run name, its
    files ('spectra_data') are the sources of the run and its database (db, db_version and taxonomy of the
    search parameters, see databaseFromLegacy()) is a database of the run, so none of them are settings.
  */
  static IdentificationData::RunSettings settingsFromLegacy(const ProteinIdentification& proteins);
  /// The database that legacy search parameters name.
  static IdentificationData::Database databaseFromLegacy(const SearchParameters& search);
  /**
    @brief A legacy protein run without proteins, identifier and files, from the settings and the first database of @p run

    If import kept a legacy protein run without an inference result although its score type was not that of the
    PSMs (e.g. empty, or the search engine score after rescoring), the settings metadata
    'identification:legacy_protein_score_type' and 'identification:legacy_protein_higher_score_better' hold it,
    and the result takes them as its score type and direction.
  */
  static ProteinIdentification settingsToLegacy(const IdentificationData::Run& run);
  /**
    @brief Keep the score type of the legacy protein run of @p run when its primary score changes

    Without an inference result, export gives the legacy protein run the score type of the primary score (as search
    engines write it), unless import recorded another one. Code that changes the primary score of a run (e.g. to
    q-values) calls this before, so the legacy protein run keeps its score type, as legacy rescoring keeps it.
  */
  static void keepLegacyProteinScoreType(IdentificationData::Run& run);
  /// The identifier of the legacy protein run of @p run, which export gives its peptide identifications
  static std::string legacyIdentifier(const IdentificationData::Run& run);
  /**
    @brief Give every match of every run of @p data a new primary score, as legacy rescoring does

    Every match gets the value of @p value (from its run, the match and its current primary score) in the score
    @p definition, which becomes the primary score. As legacy rescoring keeps the previous score as metadata of the
    hits, it becomes the meta value "<name><previous_suffix>" (or @p previous_meta, if given) of every match and is
    removed as a score (unless it is @p definition). With @p keep_different, an existing meta value of that name with a
    different value stays and the previous score becomes "<that name>~" (as IDScoreSwitcherAlgorithm::switchScores()).
    The legacy protein runs keep their score type (keepLegacyProteinScoreType()). Nothing changes for data without a
    primary score.
  */
  static void replacePrimaryScore(IdentificationData& data, const IdentificationData::ScoreDefinition& definition,
                                  const std::function<double(const IdentificationData::Run&, const IdentificationData::Match&, double)>& value,
                                  const std::string& previous_suffix, bool keep_different = false, const std::string& previous_meta = "");
  /// As above, but the previous score is not kept (e.g. when it is converted, as posterior error probabilities to posterior probabilities)
  static void replacePrimaryScore(IdentificationData& data, const IdentificationData::ScoreDefinition& definition,
                                  const std::function<double(const IdentificationData::Run&, const IdentificationData::Match&, double)>& value);
  /**
    @brief The protein hits of an inference result, completed from the database sequences of its input runs

    An import keeps a protein's sequence, description and metadata once, in the database sequences of its run;
    the hits of its inference result hold the inference values (score, rank, coverage, modifications) and the
    target/decoy state. The completed hits add what a hit does not have itself from the database sequence of its
    qualified accession, as the legacy protein run had it.
  */
  static std::vector<ProteinHit> proteinHits(const IdentificationData& data, const IdentificationData::InferenceResult& result);
  /// The protein hits of the database sequences of @p run, as export writes them in its protein run without an inference result
  static std::vector<ProteinHit> proteinHits(const IdentificationData::Run& run);
  /// The legacy protein run of @p run, as export writes it without an inference result (losses are not reported)
  static ProteinIdentification proteinRun(const IdentificationData::Run& run);
  /**
    @brief The inference result for protein inference over all peptide identifications of @p data, which export
    writes in one protein run

    This is the inference result that covers the peptide runs (e.g. of ConsensusMapMergerAlgorithm::mergeAllIDRuns()),
    with its protein hits completed (proteinHits()), or for a single peptide run without one, a new inference result
    with the legacy protein run of the run (proteinRun()) and the run as input. Its protein run has the identifier that
    export gives it. Inference stores its result with storeInferenceResult().

    @return Nothing if @p data has no peptide runs
    @throw Exception::InvalidParameter if the peptide identifications are not in one protein run
  */
  static std::optional<IdentificationData::InferenceResult> pooledInferenceResult(const IdentificationData& data);
  /// Add @p result to @p data, replacing the inference results that cover any of its input runs
  static void storeInferenceResult(IdentificationData& data, const IdentificationData::InferenceResult& result);

  /**
    @name Legacy file lists

    A legacy protein run lists its files in 'spectra_data', and its peptide identifications point into
    that list with the meta value 'id_merge_index'. Natively, the sources of a run are that list, so the
    position of a source is the legacy index and identifications need no index of their own. A source
    without a path stands for a file that is not known.
  */
  //@{
  /// Add one source per legacy file to @p run, in order, including files without identifications and repeated files.
  static void addLegacySources(IdentificationData::Run& run, const StringList& files);
  /**
    @brief The source of a legacy peptide identification in a run set up by addLegacySources() for @p n_files files

    This is the source at its 'id_merge_index', or the only file of a single-file run. Otherwise its file
    is not known, and the result is the first source without a path, which is added if there is none.

    @throw Exception::InvalidParameter if 'id_merge_index' is not an index into the @p n_files files
  */
  static IdentificationData::SourceId legacySource(IdentificationData::Run& run, Size n_files, const PeptideIdentification& item);
  /// The legacy file list of @p run: the paths of its sources that name a file, in order.
  static StringList legacyFiles(const IdentificationData::Run& run);
  //@}
};
} // namespace OpenMS
