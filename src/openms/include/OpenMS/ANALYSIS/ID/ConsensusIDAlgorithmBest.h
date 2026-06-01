// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmIdentity.h>

namespace OpenMS
{
  /**
    @brief Calculates a consensus from multiple ID runs by taking the best search score.

    @htmlinclude OpenMS_ConsensusIDAlgorithmBest.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmBest :
    public ConsensusIDAlgorithmIdentity
  {
  public:
    /// Default constructor
    ConsensusIDAlgorithmBest();

  private:
    /// Not implemented
    ConsensusIDAlgorithmBest(const ConsensusIDAlgorithmBest&);

    /// Not implemented
    ConsensusIDAlgorithmBest& operator=(const ConsensusIDAlgorithmBest&);

    /// Aggregate peptide scores into one final score (by taking the best score)
    double getAggregateScore_(std::vector<double>& scores,
                                      bool higher_better) override;
  };

} // namespace OpenMS

