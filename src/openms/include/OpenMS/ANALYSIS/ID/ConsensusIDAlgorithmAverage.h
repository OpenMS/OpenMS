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
    @brief Calculates a consensus from multiple ID runs by averaging the search scores.

    @htmlinclude OpenMS_ConsensusIDAlgorithmAverage.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmAverage :
    public ConsensusIDAlgorithmIdentity
  {
  public:
    /// Default constructor
    ConsensusIDAlgorithmAverage();

  private:
    /// Not implemented
    ConsensusIDAlgorithmAverage(const ConsensusIDAlgorithmAverage&);

    /// Not implemented
    ConsensusIDAlgorithmAverage& operator=(const ConsensusIDAlgorithmAverage&);

    /// Aggregate peptide scores into one final score (by averaging)
    double getAggregateScore_(std::vector<double>& scores,
                                      bool higher_better) override;
  };

} // namespace OpenMS

