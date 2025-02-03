// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

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
        @brief Estimates the peak type of the peaks in the iterator range based on intensity characteristics of up to five maxima
        
        We estimate profile vs. centroided by looking at highest five peaks in the spectrum.
        If at least two neighbouring sampling points to either side of a local maximum are found within 1 Th,
        the peak is considered a profile peak. The intensities need to decline on both shoulders.
        
        All sampling points successfully assigned as shoulder points are not considered for searching the next highest peak
        (except local minima at the end of a shoulder, which can be used for shoulders from left and right).
        
        If 5 peaks or 50% of total spectral intensity has been looked at, we compare the number of peaks hithero classified 
        as centroided (C) vs profile (P).
        If P / (C+P) > 0.75, the spectrum is considered profile; centroided otherwise.

        @note if there are less than 5 peaks in the iterator range SpectrumSettings::UNKNOWN is returned
       */
    template <typename PeakConstIterator>
    static SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end)
    {
      typedef typename PeakConstIterator::value_type PeakT;
      // abort if there are less than 5 peak in the iterator range
      if (end - begin < 5)
      {
        return SpectrumSettings::UNKNOWN;
      }

      const int max_peaks = 5; // maximal number of peaks we are looking at
      int profile_evidence = 0; // number of peaks found to be profile
      int centroid_evidence = 0; // number of peaks found to be centroided

      // copy data, since we need to modify
      std::vector<PeakT> data(begin, end);
      // total intensity of spectrum
      double total_int = std::accumulate(begin, end, 0.0, [](double int_, const PeakT& p) { return int_ + p.getIntensity(); } );
      double explained_int = 0;
      // get the 5 highest peaks
      for (int i = 0; i < max_peaks; ++i)
      {
        // stop if we explained +50% of all intensity 
        // (due to danger of interpreting noise - usually wrongly classified as centroided data)
        if (explained_int > 0.5 * total_int) break;
        
        double int_max = 0;
        Size idx = std::numeric_limits<Size>::max();
        // find highest peak position
        for (Size i = 0; i < data.size(); ++i)
        {
          if (data[i].getIntensity() > int_max)
          {
            int_max = data[i].getIntensity();
            idx = i;
          }
        } 
        // no more peaks
        if (idx == std::numeric_limits<Size>::max()) break;

        // check left and right peak shoulders and count number of sample points
        typedef typename std::vector<PeakT>::iterator PeakIterator; // non-const version, since we need to modify the peaks
        PeakIterator it_max = data.begin() + idx;
        PeakIterator it = it_max;
        double int_last = int_max;
        while (it != data.begin() 
               && it->getIntensity() <= int_last        // at most 100% of last sample point
               && it->getIntensity() > 0 
               && (it->getIntensity() / int_last) > 0.1 // at least 10% of last sample point
               && it->getMZ() + 1 > it_max->getMZ())    // at most 1 Th away
        {
          int_last = it->getIntensity();
          explained_int += int_last;
          it->setIntensity(0); // remove peak from future consideration
          --it;
        }
        // if the current point is rising again, restore the intensity of the 
        // previous 'sink' point (because it could belong to a neighbour peak)
        // e.g. imagine intensities: 1-2-3-2-1-2-4-2-1. We do not want to destroy the middle '1'
        if (it->getIntensity() > int_last) (it+1)->setIntensity(int_last);

        //std::cerr << "  Peak candidate: " << it_max->getMZ() << " ...";
        bool break_left = false;
        if (it_max - it < 2+1)  // 'it' does not fulfill the conditions, i.e. does not count
        { // fewer than two sampling points on left shoulder
          //std::cerr << " break left " << it_max - it << " points\n";
          break_left = true;
          // do not end loop here.. we still need to clean up the right side
        }
        it_max->setIntensity(int_max); // restore center intensity
        explained_int -= int_max;
        it = it_max;
        int_last = int_max;
        while (it != data.end() 
               && it->getIntensity() <= int_last        // at most 100% of last sample point
               && it->getIntensity() > 0 
               && (it->getIntensity() / int_last) > 0.1 // at least 10% of last sample point
               && it->getMZ() - 1 < it_max->getMZ())    // at most 1 Th away
        {
          int_last = it->getIntensity();
          explained_int += int_last;
          it->setIntensity(0); // remove peak from future consideration
          ++it;
        }
        // if the current point is rising again, restore the intensity of the 
        // previous 'sink' point (because it could belong to a neighbour peak)
        // e.g. imagine intensities: 1-2-4-2-1-2-3-2-1. We do not want to destroy the middle '1'
        // (note: the sequence is not identical to the one of the left shoulder)
        if (it != data.end() && it->getIntensity() > int_last) (it-1)->setIntensity(int_last);

        if (break_left || it - it_max < 2+1)  // 'it' does not fulfill the conditions, i.e. does not count
        { // fewer than two sampling points on right shoulder
          //std::cerr << " break right " << it - it_max << " points\n";
          ++centroid_evidence;
          continue;
        }
        // peak has at least two sampling points on either side within 1 Th
        ++profile_evidence;
        //std::cerr << " PROFILE " << it - it_max << " points right\n";

      }
      
      float evidence_ratio = profile_evidence / float(profile_evidence + centroid_evidence);
      //std::cerr << "--> Evidence ratio: " << evidence_ratio; 

      if (evidence_ratio > 0.75) // 80% are profile
      {
        //std::cerr << "  PROFILE\n";
        return SpectrumSettings::PROFILE;
      }
      else
      { 
        //std::cerr << "  CENTROID\n";
        return SpectrumSettings::CENTROID;
      }
    }
    /**
    Below code is left for reference, for things which do not work across instrument classes, resolutions and m/z ranges 
    (in case someone wants to try and improve it)

    // this code does not work reliably, mainly on nearly-empty spectra.
    // one could just use the median of minimal distances, but that requires a magical cutoff parameter, which is hard to find
    // and will fail for extreme cases, e.g. very high m/z (where resolution is decreased and thus inter-peak spacing is large for profile data)
    // Looking at quartiles of the distribution of inter-peak spacing also does not work reliable, due to variable sampling distances over the m/z range.

    Profile distances are unimodal (with very few outliers), 
    aggregating around the sampling-interval of the instrument, e.g. 0.002 Th.
    We restrict the search to the first 100 peaks, since sampling intervals can increase drastically over m/z, thus
    distorting the unimodal model, e.g. for Orbitrap 0.0006@200 Th to 0.1@4000 Th.

    Centroided distances are multimodal, aggregating at 1, 1/2, 1/3 etc and surprisingly often near 0 as well.

    Comparing the distance between the first and third quantile of the distance distribution gives an indication
    on the type of data. The Q1 vs. Q3 only differs by a few percent [(Q1-Q3)/Q1*100 ~ 4%] we are surely looking at profile data.
    On the other hand, observing large differences (~ factor 10-100) indicates centroided data.
    We set the threshold to factor of 5, i.e. if (Q3-Q1)/Q1 < 5, it's profile data.

    template <typename PeakConstIterator>
    static SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end)
    {
      const int MAX_SAMPLED_DISTANCES = 100;
      const double MAX_MZ_WINDOW = 300; // inspect only a small window, otherwise sampling distance between raw peaks will change too much
                                        // abort if there are less than 5 peak in the iterator range
      if (end - begin < 5)
      {
        return SpectrumSettings::UNKNOWN;
      }

      int count(0);

      std::vector<double> distances;

      PeakConstIterator peak(begin);

      while (peak->getIntensity() <= 0 && peak != end)  
      { // 1st positive intensity
        ++peak;
      }

      if (peak == end)
      { // only zeros
        return SpectrumSettings::UNKNOWN;
      }

      double last_mz = peak->getMZ();
      double last_dist = std::numeric_limits<double>::max();

      for (++peak; peak != end
        && count < MAX_SAMPLED_DISTANCES
        && peak->getMZ() - begin->getMZ() < MAX_MZ_WINDOW
        ; ++peak)
      {
        if (peak->getIntensity() <= 0) continue;

        double dist = peak->getMZ() - last_mz;
        distances.push_back(std::min(last_dist, dist)); // min distances to either side
        ++count;
        last_mz = peak->getMZ();
        last_dist = dist;
      }

      if (count < 4) // at least 4 distances for non-zero(!) intensity peaks
      {
        if (peak != end) return estimateType(peak, end); // try further to the right
        else return SpectrumSettings::UNKNOWN;
      }

      double q1 = Math::quantile1st(distances.begin(), distances.end(), false);
      double q3 = Math::quantile3rd(distances.begin(), distances.end(), true);
      for (const auto& i : distances) std::cerr << i << ";";
      std::cerr <<"\t";

      if ((q3-q1) < q1*5) // q1 and q3 are roughly equal
      {
        return SpectrumSettings::PROFILE;
      }
      else
      {
        return SpectrumSettings::CENTROID;
      }
    }*/

  };

} // namespace OpenMS

