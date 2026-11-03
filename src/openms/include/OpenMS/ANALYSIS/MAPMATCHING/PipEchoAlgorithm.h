// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/OpenMSConfig.h>
#include <vector>

namespace OpenMS
{

  /**
    @brief Match-between-runs feature grouping via the PIP-ECHO algorithm.

    Groups features that represent the same identified peptide across multiple
    retention-time-aligned LC-MS runs of a single fraction. Identifications are
    transferred to runs in which a feature was detected but not identified
    (match-between-runs); candidate transfers are scored and filtered at a
    controlled false-discovery rate using a target/decoy model.

    A consensus feature that received transfers lists the map indices of the runs they went to in
    the meta value Constants::UserParam::MBR_TRANSFER_MAP_INDICES ("mbr_transfer_map_indices") and
    their q-values, in the same order, in Constants::UserParam::MBR_TRANSFER_QVALUES
    ("mbr_transfer_qvalues").

    @htmlinclude OpenMS_PipEchoAlgorithm.parameters

    @ingroup FeatureGrouping
  */
class OPENMS_DLLAPI PipEchoAlgorithm : public FeatureGroupingAlgorithm
{
public:
  /// Constructor.
  PipEchoAlgorithm();

  PipEchoAlgorithm(PipEchoAlgorithm&&) = delete;

  PipEchoAlgorithm& operator=(PipEchoAlgorithm&&) = delete;

  /// Copy constructor.
  PipEchoAlgorithm(const PipEchoAlgorithm& other) = delete;

  /// Assignment operator.
  PipEchoAlgorithm& operator=(const PipEchoAlgorithm& source) = delete;

  /// Destructor.
  ~PipEchoAlgorithm() override;

  /**
    @brief Groups corresponding features across the feature maps of one fraction.

    Features representing the same identified peptide are grouped together and,
    where possible, identifications are transferred between runs.

    @param features The input feature maps (same fraction) to link.
    @param consensus The resulting consensus map of grouped features.

    @pre The input maps must be retention-time aligned (e.g. via a map-alignment
         algorithm) and originate from the same fraction.
  */
  void group(const std::vector<FeatureMap>& features, ConsensusMap& consensus) override;
};

} // namespace OpenMS
