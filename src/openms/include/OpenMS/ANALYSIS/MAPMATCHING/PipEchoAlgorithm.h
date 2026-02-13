// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
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
  ~PipEchoAlgorithm();

  /**
   * Group together features from multiple LC-MS runs that represent
   * the same identified peptides.
   *
   * The given feature maps must be aligned prior to using this
   * function and come from the same fraction.
   */
  void group(const std::vector<FeatureMap>& features, ConsensusMap& consensus) override;
};

} // namespace OpenMS
