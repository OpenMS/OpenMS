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
    @brief Calculates a consensus from multiple ID runs based on the ranks of the search hits.

    @htmlinclude OpenMS_ConsensusIDAlgorithmRanks.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmRanks :
    public ConsensusIDAlgorithmIdentity
  {
  public:
    /// Default constructor
    ConsensusIDAlgorithmRanks();

  private:
    /// Number of ID runs for current analysis
    Size current_number_of_runs_;

    /// Number of considered hits for current analysis
    Size current_considered_hits_;

    /// Not implemented
    ConsensusIDAlgorithmRanks(const ConsensusIDAlgorithmRanks&);

    /// Not implemented
    ConsensusIDAlgorithmRanks& operator=(const ConsensusIDAlgorithmRanks&);

    /// Assign peptide scores based on search ranks
    void preprocess_(std::vector<PeptideIdentification>& ids) override;

    /// Aggregate peptide scores into one final score (by averaging ranks)
    double getAggregateScore_(std::vector<double>& scores,
                                      bool higher_better) override;
   };

} // namespace OpenMS

