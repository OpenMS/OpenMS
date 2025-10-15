// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS {
  class OPENMS_DLLAPI PipEcho
    : public FeatureGroupingAlgorithm
  {
 public:

    /// Constructor.
    PipEcho();

    PipEcho(PipEcho&&) = delete;

    PipEcho& operator=(PipEcho&&) = delete;

    /// Copy constructor.
    PipEcho(const PipEcho& other) = delete;

    /// Assignment operator.
    PipEcho& operator=(const PipEcho& source) = delete;

    /// Destructor.
    ~PipEcho();

    /**
     * Group together features from multiple LC-MS runs that represent
     * the same identified peptides.
     *
     * The given feature maps must be aligned prior to using this
     * function and come from the same fraction.
     */
    void group(
      const std::vector<FeatureMap>& features,
      ConsensusMap& consensus
    ) override;
  };
}
