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
#include <optional>
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
    @brief Move the peptide identifications of a map into its identification data, for code that edits identification data

    @return Whether @p map had peptide identifications, i.e. whether to move them back afterwards (exportFeatureIDs())
    @throw Exception::InvalidParameter if @p map has peptide identifications and identification data
  */
  static bool moveToIdentificationData(FeatureMap& map);
  static bool moveToIdentificationData(ConsensusMap& map);
  ///@}
};
}
