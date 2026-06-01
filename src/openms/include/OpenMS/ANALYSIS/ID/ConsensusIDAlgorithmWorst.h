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
    @brief Calculates a consensus from multiple ID runs by taking the worst search score (conservative approach).

    @htmlinclude OpenMS_ConsensusIDAlgorithmWorst.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmWorst :
    public ConsensusIDAlgorithmIdentity
  {
  public:
    /// Default constructor
    ConsensusIDAlgorithmWorst();

  private:
    /// Not implemented
    ConsensusIDAlgorithmWorst(const ConsensusIDAlgorithmWorst&);

    /// Not implemented
    ConsensusIDAlgorithmWorst& operator=(const ConsensusIDAlgorithmWorst&);

    /// Aggregate peptide scores into one final score (by taking the worst score)
    double getAggregateScore_(std::vector<double>& scores,
                                      bool higher_better) override;
  };

} // namespace OpenMS

