// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>
#include <iomanip>

namespace OpenMS
{
  float Qscore::getQscore(const PeakGroup* pg)
  {
    if (pg->empty())
    { // all zero
      return .0f;
    }
    // the weights for per cosine, SNR, PPM error, charge score, and intercept.
    // Cos    -8.9494
    // SNR    -2.2977
    // PPM error           0.3269
    // charge score        -3.1663
    // Intercept    12.3131
    // const std::vector<double> weights({-8.9494, -2.2977, 0.3269, -3.1663, 12.3131});
    const std::vector<double> weights({-2.2833, -3.2881, 0, 0, 4.5425});
    double score = weights.back();
    auto fv = toFeatureVector_(pg);

    for (Size i = 0; i < weights.size() - 1; i++)
    {
      score += fv[i] * weights[i];
    }
    float qscore = 1.0f / (1.0f + (float)exp(score));

    return qscore;
  }

  std::vector<double> Qscore::toFeatureVector_(const PeakGroup* pg)
  {
    std::vector<double> fvector(4); // length of weights vector - 1, excluding the intercept weight.

    double a = pg->getIsotopeCosine();
    double d = 1;
    int index = 0;
    fvector[index++] = (log2(a + d));

    a = pg->getSNR();
    fvector[index++] = (log2(d + a / (d + a)));

    a = pg->getAvgPPMError();
    fvector[index++] = (log2(d + a / (d + a)));

    a = pg->getChargeScore();
    fvector[index++] = (log2(a + d));

    return fvector;
  }
} // namespace OpenMS
