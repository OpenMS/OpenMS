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
#include <OpenMS/METADATA/ID/MapPeptideIdentifications.h>
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
  /// The identifications of a map as peptide and protein identifications, apart from the map (see exportMapIdentifications())
  using MapPeptideIdentifications = OpenMS::MapPeptideIdentifications;
  /**
    @brief The identifications of @p map (its identification data and the links of its features) as peptide and
    protein identifications, as a writer of a format holding those writes them
  */
  static MapPeptideIdentifications exportMapIdentifications(const FeatureMap& map);
  static MapPeptideIdentifications exportMapIdentifications(const ConsensusMap& map);
  /**
    @brief Set the identification data of @p map and the links of its features (and subordinates) from peptide and
    protein identifications (e.g. those of a file), replacing the identifications of @p map

    @throw Exception::InvalidParameter if a path refers to no feature of @p map, or for identifications that
    identification data cannot hold (see IdentificationDataAdapter::importLegacyMap())
  */
  static void importMapIdentifications(FeatureMap& map, const MapPeptideIdentifications& identifications);
  static void importMapIdentifications(ConsensusMap& map, const MapPeptideIdentifications& identifications);

  /// @name Maps for code that works on identification data
  ///@{
  /**
    @brief The protein identifications of a map, as export writes them (exportMapIdentifications()), e.g. the protein
    run of an inference result
  */
  static std::vector<ProteinIdentification> proteinIdentifications(const FeatureMap& map);
  static std::vector<ProteinIdentification> proteinIdentifications(const ConsensusMap& map);
  /**
    @brief Append the features of @p source to @p target with the identifications they link, e.g. the parts of a map
    that parallel producers built separately

    A run of @p source is appended to the run of @p target with the same identifier (see IdentificationData::Run::append()),
    or added if @p target has none; the links of the appended features (and their subordinates) follow. Identifications
    of @p source that no feature links are appended as well. Inference results of @p source are not.

    @throw Exception::InvalidValue if runs with the same identifier differ in molecule kind or score definitions
  */
  static void appendFeatures(FeatureMap& target, const FeatureMap& source);
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
    @brief The maps with each identification run in one map only (to combine their identifications, e.g. when
    grouping features)

    A run that is also in an earlier map gets a new UUID (see makeRunsDistinct()). Maps are copied only if one needs
    to change.

    @return @p maps, or @p converted holding all maps
  */
  static const std::vector<FeatureMap>& withDistinctRuns(const std::vector<FeatureMap>& maps, std::vector<FeatureMap>& converted);
  static const std::vector<ConsensusMap>& withDistinctRuns(const std::vector<ConsensusMap>& maps, std::vector<ConsensusMap>& converted);
  /**
    @brief The identifications of @p map as peptide identifications (see exportMapIdentifications()) whose peptide hits
    name their match (see matchReference()), for algorithms that work on peptide hits, e.g. protein inference
  */
  static MapPeptideIdentifications exportWithMatchReferences(const ConsensusMap& map);
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
