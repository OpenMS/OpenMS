// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
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

        @note if there are fewer than 5 peaks in the iterator range SpectrumSettings::UNKNOWN is returned
       */
    template <typename PeakConstIterator>
    SpectrumSettings::SpectrumType estimateType(const PeakConstIterator & begin, const PeakConstIterator & end) const
    {
      const Size MAX_SAMPLED_DISTANCES = 1000;
      const double DISTANCE_VARIANCE_THRESHOLD = 0.5;

      // abort if there are less than 5 peak in the iterator range
      if (end - begin < 5)
      {
        return SpectrumSettings::UNKNOWN;
      }

      double count(0);

      std::vector<double> distances;

      PeakConstIterator peak(begin);

      for (; peak->getIntensity() <= 0 && peak != end - 2; ++peak)  // 1st positive intensity
      {
      }

      double scnd_last_mz(peak->getMZ());

      for (++peak; peak->getIntensity() <= 0 && peak != end - 1; ++peak) // 2nd positive intensity
      {
      }

      double last_mz(peak->getMZ());

      double last_dist(last_mz - scnd_last_mz);

      for (++peak; peak != end && count < MAX_SAMPLED_DISTANCES; ++peak) // max  positive intensity
      {
        if (peak->getIntensity() > 0)
        {
          double mz(peak->getMZ());
          double dist(mz - last_mz);
          distances.push_back(std::min(last_dist, dist)); // min distances
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

      double mean(std::accumulate(distances.begin(), distances.end(), 0) / count); // sum/size

      // calculate variance
      double variance(0);
      for (std::vector<double>::iterator value = distances.begin(); value != distances.end(); ++value)
      {
        double delta = (*value - mean);
        variance += delta * delta;
      }
      variance /= count - 1;

      // calculate stdev
      double standard_deviation(std::sqrt(variance));

      if (standard_deviation < DISTANCE_VARIANCE_THRESHOLD)
      {
        return SpectrumSettings::PROFILE;
      }
      else
      {
        return SpectrumSettings::CENTROID;
      }
    }

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEAKTYPEESTIMATOR_H
