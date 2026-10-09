// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <functional>
#include <optional>
#include <set>
#include <string>
#include <vector>
namespace OpenMS
{
class FeatureMap;
class ConsensusMap;
/** @brief Explicit adapters for owning identifications, FASTA databases and map annotations. */
class OPENMS_DLLAPI IdentificationDataConverter
{
public:
  static void importIDs(IdentificationData& data, const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides);
  static void exportIDs(const IdentificationData& data,
                        std::vector<ProteinIdentification>& proteins,
                        PeptideIdentificationList& peptides,
                        bool export_ids_wo_scores = false);
  static MzTab exportMzTab(const IdentificationData& data);
  /// Add @p database to @p run (or find an equal one) and its FASTA entries as database sequences; returns the database.
  static IdentificationData::DatabaseId importSequences(IdentificationData::Run& run, const IdentificationData::Database& database,
                                                        const std::vector<FASTAFile::FASTAEntry>& fasta, const std::string& decoy_pattern = "");
  static void exportSequenceEvidence(const std::vector<IdentificationData::SequenceEvidence>& evidence, PeptideHit& hit);
  static void importFeatureIDs(FeatureMap& features, bool clear_original = true);
  static void exportFeatureIDs(FeatureMap& features, bool clear_original = true);
  static void importConsensusIDs(ConsensusMap& consensus, bool clear_original = true);
  static void exportConsensusIDs(ConsensusMap& consensus, bool clear_original = true);

  /// @name Maps for code that works on identification data
  ///@{
  /// Whether @p map has peptide or protein identifications (of features, subordinates or unassigned)
  static bool hasPeptideIdentifications(const FeatureMap& map);
  static bool hasPeptideIdentifications(const ConsensusMap& map);
  /**
    @brief Read access to the identifications of a map as identification data

    @return @p map if it has no peptide identifications, else @p converted: a copy of @p map with its peptide
    identifications moved into its identification data (importFeatureIDs())

    @throw Exception::InvalidParameter if @p map has peptide identifications and identification data
  */
  static const FeatureMap& withIdentificationData(const FeatureMap& map, std::optional<FeatureMap>& converted);
  static const ConsensusMap& withIdentificationData(const ConsensusMap& map, std::optional<ConsensusMap>& converted);
  /**
    @brief Read access to the identifications of a map as peptide identifications, for code that works on them (e.g.
    writers of formats that hold peptide identifications)

    The converse of withIdentificationData().

    @return @p map if it has no identification data or has peptide identifications, else @p exported: a copy of @p map
    with its identification data exported as peptide identifications (exportFeatureIDs())
  */
  static const FeatureMap& withPeptideIdentifications(const FeatureMap& map, std::optional<FeatureMap>& exported);
  static const ConsensusMap& withPeptideIdentifications(const ConsensusMap& map, std::optional<ConsensusMap>& exported);
  /**
    @brief The protein identifications of a map: its own, or for a map with identification data (and none of its own)
    those that export writes (exportFeatureIDs()), e.g. the protein run of an inference result
  */
  static std::vector<ProteinIdentification> proteinIdentifications(const FeatureMap& map);
  static std::vector<ProteinIdentification> proteinIdentifications(const ConsensusMap& map);
  /**
    @brief Move the peptide identifications of a map into its identification data, for code that edits identification data

    @return Whether @p map had peptide identifications, i.e. whether to move them back afterwards (exportFeatureIDs())
    @throw Exception::InvalidParameter if @p map has peptide identifications and identification data
  */
  static bool moveToIdentificationData(FeatureMap& map);
  static bool moveToIdentificationData(ConsensusMap& map);
  /**
    @brief Edit the identifications of a map as identification data

    Runs @p edit on @p map. A map with peptide identifications has them moved into its identification data for the
    edit (moveToIdentificationData()) and back afterwards (exportFeatureIDs()/exportConsensusIDs()), also if @p edit
    throws. Its protein identification runs keep their objects if their number does not change, so references to
    them (e.g. to the protein run of an inference result) stay valid.

    @throw Exception::InvalidParameter if @p map has peptide identifications and identification data
  */
  static void editAsIdentificationData(FeatureMap& map, const std::function<void(FeatureMap&)>& edit);
  static void editAsIdentificationData(ConsensusMap& map, const std::function<void(ConsensusMap&)>& edit);
  /**
    @brief Give the identification runs of @p map that are in @p taken new UUIDs, so that its identifications stay apart
    from those of the maps that took them when they are combined (e.g. when grouping features)

    The links of the features (and subordinates) and the inputs of inference results follow the new UUIDs. The UUIDs of
    the runs of @p map are added to @p taken.

    @return Whether a run got a new UUID
  */
  static bool makeRunsDistinct(FeatureMap& map, std::set<std::string>& taken);
  static bool makeRunsDistinct(ConsensusMap& map, std::set<std::string>& taken);
  /**
    @brief The maps with their identifications as identification data, each identification run in one map only (to
    combine their identifications, e.g. when grouping features)

    Maps with peptide identifications are converted (see withIdentificationData()); a run that is also in an earlier
    map gets a new UUID (see makeRunsDistinct()). Maps are copied only if one needs to change.

    @return @p maps, or @p converted holding all maps
    @throw Exception::InvalidParameter if a map has peptide identifications and identification data
  */
  static const std::vector<FeatureMap>& withIdentificationData(const std::vector<FeatureMap>& maps, std::vector<FeatureMap>& converted);
  static const std::vector<ConsensusMap>& withIdentificationData(const std::vector<ConsensusMap>& maps, std::vector<ConsensusMap>& converted);
  /**
    @brief A copy of @p map with its identifications as peptide identifications (exportConsensusIDs()) whose peptide
    hits name their match (see matchReference()), for algorithms that work on peptide hits, e.g. protein inference

    @throw Exception::InvalidParameter if @p map has peptide identifications
  */
  static ConsensusMap exportWithMatchReferences(const ConsensusMap& map);
  /// The match that a peptide hit of exportWithMatchReferences() stands for, if it names one
  static std::optional<IdentificationData::MatchReference> matchReference(const PeptideHit& hit);
  /**
    @brief Take what an algorithm changed in the peptide hits of exportWithMatchReferences() back to their matches in
    @p data

    The sequence evidence of a match keeps the proteins that its hits still refer to. With @p scores, the primary
    score of a match becomes the score of its (first) hit. Matches without hits in @p peptides stay as they are.
  */
  static void updateReferencedMatches(IdentificationData& data, const PeptideIdentificationList& peptides, bool scores);
  ///@}
};
}
