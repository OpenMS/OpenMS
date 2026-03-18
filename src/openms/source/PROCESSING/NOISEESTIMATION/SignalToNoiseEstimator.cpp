// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <random>

using namespace std;

namespace OpenMS
{

  float estimateNoiseFromRandomScans(const MSExperiment& exp, const UInt ms_level, const UInt n_scans, const double percentile)
  {
    vector<Size> spec_indices;
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (exp[i].getMSLevel() == ms_level && !exp[i].empty())
      {
        spec_indices.push_back(i);
      }
    }

    if (spec_indices.empty()) return 0.0f;

    std::default_random_engine generator(time(nullptr));
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    float noise = 0.0;
    UInt count = 0;
    vector<float> tmp;
    while (count++ < n_scans)
    {
      UInt scan = (UInt)(distribution(generator) * (spec_indices.size() - 1));
      tmp.clear();
      for (const auto& peak : exp[scan])
      {
        tmp.push_back(peak.getIntensity());
      }
      Size idx = tmp.size() * percentile / 100.0;
      std::nth_element(tmp.begin(), tmp.begin() + idx, tmp.end());
      noise += tmp[idx];
    }
    return noise / n_scans;
  }

}
