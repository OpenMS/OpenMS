// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmWorst.h>

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmWorst::ConsensusIDAlgorithmWorst()
  {
    setName("ConsensusIDAlgorithmWorst"); // DefaultParamHandler
  }


  double ConsensusIDAlgorithmWorst::getAggregateScore_(vector<double>& scores,
                                                       bool higher_better)
  {
    if (higher_better)
    {
      return *min_element(scores.begin(), scores.end());
    }
    return *max_element(scores.begin(), scores.end());
  }

} // namespace OpenMS
