// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmBest.h>

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmBest::ConsensusIDAlgorithmBest()
  {
    setName("ConsensusIDAlgorithmBest"); // DefaultParamHandler
  }


  double ConsensusIDAlgorithmBest::getAggregateScore_(vector<double>& scores,
                                                      bool higher_better)
  {
    if (higher_better)
    {
      return *max_element(scores.begin(), scores.end());
    }
    return *min_element(scores.begin(), scores.end());
  }

} // namespace OpenMS
