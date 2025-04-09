// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski  $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>

#include <algorithm>
#include <numeric>
#include <functional>
#include <stdexcept>

namespace OpenSwath
{
  void normalize(
    const std::vector<double> & intensities,
    double normalizer,
    std::vector<double> & normalized_intensities)
  {
    //compute total intensities
    //normalize intensities
    normalized_intensities.resize(intensities.size());
    if (normalizer > 0)
    {
      std::transform(intensities.begin(), intensities.end(), normalized_intensities.begin(),
                     [&normalizer](double val)
                     {
                      return val / normalizer;
                     });
    }
  }


  double dotprodScoring(std::vector<double> intExp, std::vector<double> theorint)
  {
    for (unsigned int i = 0; i < intExp.size(); ++i)
    {
      intExp[i] = sqrt(intExp[i]);
      theorint[i] = sqrt(theorint[i]);
    }

    double intExptotal = norm(intExp.begin(), intExp.end());
    double intTheorTotal = norm(theorint.begin(), theorint.end());
    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::dotProd(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

  double manhattanScoring(std::vector<double> intExp, std::vector<double> theorint)
  {

    for (unsigned int i = 0; i < intExp.size(); ++i)
    {
      intExp[i] = sqrt(intExp[i]);
      theorint[i] = sqrt(theorint[i]);
      //std::transform(intExp.begin(), intExp.end(), intExp.begin(), sqrt);
      //std::transform(theorint.begin(), theorint.end(), theorint.begin(), sqrt);
    }

    double intExptotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
    double intTheorTotal = std::accumulate(theorint.begin(), theorint.end(), 0.0);
    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

}
