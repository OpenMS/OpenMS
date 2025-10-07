// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
  class OPENMS_DLLAPI IonIdentityMolecularNetworking
  {
    public:
      /// Annotate ConsensusMap for ion identity molecular networking (IIMN) workflow by GNPS.
      /// Adds meta values Constants::UserParams::IIMN_ROW_ID (unique index for each feature), Constants::UserParams::IIMN_ADDUCT_PARTNERS (related features row IDs)
      /// and Constants::UserParams::IIMN_ANNOTATION_NETWORK_NUMBER (all related features with different adduct states) get the same network number).
      /// This method requires the features annotated with the Constants::UserParams::IIMN_LINKED_GROUPS meta value.
      ///  If at least one of the features has an annotation for Constants::UserParam::IIMN_LINKED_GROUPS, annotate ConsensusMap for IIMN.
      static void annotateConsensusMap(ConsensusMap& consensus_map);

      /// Write supplementary pair table (csv file) from a consensusXML file with edge annotations for connected features. Required for GNPS IIMN.
      /// The table contains the columns "ID 1" (row ID of first feature), "ID 2" (row ID of second feature), "EdgeType" (MS1/2 annotation),
      /// "Score" (the number of direct partners from both connected features) and "Annotation" (adducts and delta m/z between two connected features).
      static void writeSupplementaryPairTable(const ConsensusMap& consensus_map, const String& output_file);
  };
} // closing namespace OpenMS