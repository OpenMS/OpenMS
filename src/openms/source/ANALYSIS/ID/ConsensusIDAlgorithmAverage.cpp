// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmAverage.h>

#include <numeric> // for "accumulate"

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmAverage::ConsensusIDAlgorithmAverage()
  {
    setName("ConsensusIDAlgorithmAverage"); // DefaultParamHandler
  }


  double ConsensusIDAlgorithmAverage::getAggregateScore_(
    vector<double>& scores, bool /* higher_better */)
  {
    double sum_scores = accumulate(scores.begin(), scores.end(), 0.0);
    return sum_scores / scores.size();
  }

} // namespace OpenMS
