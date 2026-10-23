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

#include <algorithm>
#include <random>

using namespace std;

namespace OpenMS
{

  float estimateNoiseFromRandomScans(const MSExperiment& exp, const UInt ms_level, const UInt n_scans, const double percentile)
  {
    if (n_scans == 0)
    {
      return 0.0f;
    }

    vector<Size> spec_indices;
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (exp[i].getMSLevel() == ms_level && !exp[i].empty())
      {
        spec_indices.push_back(i);
      }
    }

    if (spec_indices.empty())
    {
      return 0.0f;
    }

    std::default_random_engine generator(time(nullptr));
    std::uniform_int_distribution<Size> distribution(0, spec_indices.size() - 1);
    const double bounded_percentile = std::clamp(percentile, 0.0, 100.0);

    float noise = 0.0f;
    vector<float> tmp;
    for (UInt count = 0; count < n_scans; ++count)
    {
      // Sample only non-empty spectra of the requested MS level.
      const Size scan = spec_indices[distribution(generator)];
      tmp.clear();
      for (const auto& peak : exp[scan])
      {
        tmp.push_back(peak.getIntensity());
      }
      const Size idx = std::min(static_cast<Size>(tmp.size() * bounded_percentile / 100.0), tmp.size() - 1);
      std::nth_element(tmp.begin(), tmp.begin() + idx, tmp.end());
      noise += tmp[idx];
    }
    return noise / n_scans;
  }

}
