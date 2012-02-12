// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm, Marcel Schilling $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEAKTYPEESTIMATOR_H
#define OPENMS_FORMAT_PEAKTYPEESTIMATOR_H

#include <OpenMS/METADATA/SpectrumSettings.h>

#include <cmath>
#include <numeric>

namespace OpenMS
{
/**
  @brief Estimates if the data of a spectrum is raw data or peak data
  
   @ingroup Format
 */
class OPENMS_DLLAPI PeakTypeEstimator
{
public:
  /**
      @brief Estimates the peak type of the peaks in the iterator range based on the variance of inter-peak distances

      @note if there are fewer than 5 peaks in the iterator range SpectrumSettings::UNKOWN is returned
     */
  template <typename PeakConstIterator>
  SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end) const
  {
    const Size MAX_SAMPLED_DISTANCES = 1000;
    const DoubleReal DISTANCE_VARIANCE_THRESHOLD = 0.5;

    // abort if there are less than 5 peak in the iterator range
    if (end - begin < 5)
    {
      return SpectrumSettings::UNKNOWN;
    }

    DoubleReal count(0);

    std::vector<DoubleReal> distances;

    PeakConstIterator peak(begin);

    for(;peak->getIntensity() <= 0 && peak != end-2; ++peak)        // 1st positive intensity
    {
    }

    DoubleReal scnd_last_mz(peak->getMZ());

    for(++peak;peak->getIntensity() <= 0 && peak != end-1; ++peak)  // 2nd positive intensity
    {
    }

    DoubleReal last_mz(peak->getMZ());

    DoubleReal last_dist(last_mz - scnd_last_mz);

    for(++peak; peak != end && count < MAX_SAMPLED_DISTANCES; ++peak)  // max  positive intensity
    {
      if(peak->getIntensity() > 0)
      {
        DoubleReal mz(peak->getMZ());
        DoubleReal dist(mz - last_mz);
        distances.push_back(std::min(last_dist, dist));  // min distances
        ++count;
        scnd_last_mz = last_mz;
        last_mz = mz;
        last_dist = dist;
      }
    }

    if (count < 4) // at least 4 distances for non-zero(!) intensity peaks
    {
      return SpectrumSettings::UNKNOWN;
    }

    DoubleReal mean( std::accumulate(distances.begin(), distances.end(), 0) / count); // sum/size

    // calculate variance
    DoubleReal variance(0);
    for (std::vector<DoubleReal>::iterator value = distances.begin(); value != distances.end(); ++value)
    {
      DoubleReal delta = (*value - mean);
      variance += delta * delta;
    }
    variance /= count-1;

    // calculate stdev
    DoubleReal standard_deviation(std::sqrt(variance));

    if (standard_deviation < DISTANCE_VARIANCE_THRESHOLD)
    {
      return SpectrumSettings::RAWDATA;
    }
    else
    {
      return SpectrumSettings::PEAKS;
    }
  }

};

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEAKTYPEESTIMATOR_H
