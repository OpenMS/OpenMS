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
#include <set>
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
  /**
    @brief Export to mzTab with caller-selected optional metadata columns

    Oligonucleotide runs are exported as nucleic acid, oligonucleotide and OSM rows. Each OSM row gets an
    'opt_global_<key>' column for every key of @p observation_meta_keys and @p match_meta_keys (empty if neither the
    query nor the match has the value; the match value takes precedence), and 'opt_global_adduct' for a match with
    an adduct. Peptide runs are exported through the established classes (MzTab::exportIdentificationsToMzTab()),
    with all PSM metadata.

    @param[in] data Identification data to export
    @param[in] observation_meta_keys Metadata of queries (observations) exported as opt_global_<key>
    @param[in] match_meta_keys Metadata of matches exported as opt_global_<key>
  */
  static MzTab exportMzTab(const IdentificationData& data, const std::set<std::string>& observation_meta_keys = {},
                           const std::set<std::string>& match_meta_keys = {});
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
