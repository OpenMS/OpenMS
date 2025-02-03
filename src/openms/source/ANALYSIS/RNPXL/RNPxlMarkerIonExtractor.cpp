// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

using namespace std;

namespace OpenMS
{


RNPxlMarkerIonExtractor::MarkerIonsType RNPxlMarkerIonExtractor::extractMarkerIons(const PeakSpectrum& s, const double marker_tolerance)
{
  MarkerIonsType marker_ions;
  marker_ions["A"].push_back(make_pair(136.06231, 0.0));
  marker_ions["A"].push_back(make_pair(330.06033, 0.0));
  marker_ions["C"].push_back(make_pair(112.05108, 0.0));
  marker_ions["C"].push_back(make_pair(306.04910, 0.0));
  marker_ions["G"].push_back(make_pair(152.05723, 0.0));
  marker_ions["G"].push_back(make_pair(346.05525, 0.0));
  marker_ions["U"].push_back(make_pair(113.03509, 0.0));
  marker_ions["U"].push_back(make_pair(307.03311, 0.0));

  PeakSpectrum spec(s);
  Normalizer normalizer;
  normalizer.filterSpectrum(spec);
  spec.sortByPosition();

  // for each nucleotide with marker ions
  for (MarkerIonsType::iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
  {
    // for each marker ion of the current nucleotide
    for (Size i = 0; i != it->second.size(); ++i)
    {
      double mz = it->second[i].first;
      double max_intensity = 0;
      for (PeakSpectrum::ConstIterator sit = spec.begin(); sit != spec.end(); ++sit)
      {
        if (sit->getMZ() + marker_tolerance < mz)
        {
          continue;
        }
        if (mz < sit->getMZ() - marker_tolerance)
        {
          break;
        }
        if (fabs(mz - sit->getMZ()) < marker_tolerance)
        {
          if (max_intensity < sit->getIntensity())
          {
            max_intensity = sit->getIntensity();
          }
        }
      }
      it->second[i].second = max_intensity;
    }
  }
  return marker_ions;
}
}

