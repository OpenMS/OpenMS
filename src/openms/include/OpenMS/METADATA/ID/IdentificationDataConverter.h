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
namespace OpenMS
{
class FeatureMap;
class ConsensusMap;
/** @brief Explicit adapters for owning identifications, FASTA parents and map annotations. */
class OPENMS_DLLAPI IdentificationDataConverter
{
public:
  static void importIDs(IdentificationData& data, const std::vector<ProteinIdentification>& proteins, const PeptideIdentificationList& peptides);
  static void exportIDs(const IdentificationData& data,
                        std::vector<ProteinIdentification>& proteins,
                        PeptideIdentificationList& peptides,
                        bool export_ids_wo_scores = false);
  static MzTab exportMzTab(const IdentificationData& data);
  static void importSequences(IdentificationData::Run& run, const std::vector<FASTAFile::FASTAEntry>& fasta, const std::string& decoy_pattern = "");
  static void exportParentMatches(const std::vector<IdentificationData::ParentEvidence>& evidence, PeptideHit& hit);
  static void importFeatureIDs(FeatureMap& features, bool clear_original = true);
  static void exportFeatureIDs(FeatureMap& features, bool clear_original = true);
  static void importConsensusIDs(ConsensusMap& consensus, bool clear_original = true);
  static void exportConsensusIDs(ConsensusMap& consensus, bool clear_original = true);
};
}
