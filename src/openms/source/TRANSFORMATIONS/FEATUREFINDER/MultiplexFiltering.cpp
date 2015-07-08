// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultPeak.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <QDir>

#include <QDir>

using namespace std;
using namespace boost::math;

namespace OpenMS
{

  MultiplexFiltering::MultiplexFiltering(const MSExperiment<Peak1D>& exp_picked, const std::vector<MultiplexPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling) :
    exp_picked_(exp_picked), patterns_(patterns), peaks_per_peptide_min_(peaks_per_peptide_min), peaks_per_peptide_max_(peaks_per_peptide_max), missing_peaks_(missing_peaks), intensity_cutoff_(intensity_cutoff), mz_tolerance_(mz_tolerance), mz_tolerance_unit_(mz_tolerance_unit), peptide_similarity_(peptide_similarity), averagine_similarity_(averagine_similarity), averagine_similarity_scaling_(averagine_similarity_scaling)
  {
  }

  int MultiplexFiltering::positionsAndBlacklistFilter(MultiplexPeakPattern pattern, int spectrum, vector<double> peak_position, int peak, vector<double>& mz_shifts_actual, vector<int>& mz_shifts_actual_indices) const
  {
    // Try to find peaks at the expected m/z positions
    // loop over expected m/z shifts of a peak pattern
    for (unsigned mz_position = 0; mz_position < pattern.getMZShiftCount(); ++mz_position)
    {
      double scaling = 1;
      if (mz_position % (peaks_per_peptide_max_ + 1) == 0)
      {
        // Let us be more lenient when looking for zeroths peaks
        // i.e. allow for an increased deviation between expected m/z position and the actual one
        scaling = 2;
      }

      int index = getPeakIndex(peak_position, peak, peak_position[peak] + pattern.getMZShiftAt(mz_position), scaling);
      if (index != -1)
      {
        mz_shifts_actual.push_back(peak_position[index] - peak_position[peak]);
        mz_shifts_actual_indices.push_back(index);
      }
      else
      {
        mz_shifts_actual.push_back(std::numeric_limits<double>::quiet_NaN());
        mz_shifts_actual_indices.push_back(-1);
      }

    }

    // remove peaks which run into the next peptide
    // i.e. the isotopic peak of one peptide lies to the right of the mono-isotopic peak of the next one
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount() - 1; ++peptide)
    {
      double mz_shift_next_peptide = mz_shifts_actual[(peptide + 1) * (peaks_per_peptide_max_ + 1) + 1]; // m/z shift of the mono-isotopic peak of the following peptide
      if (!(boost::math::isnan)(mz_shift_next_peptide))
      {
        for (int isotope = 0; isotope < peaks_per_peptide_max_; ++isotope)
        {
          int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1; // index in m/z shift list
          if (mz_shifts_actual[mz_position] >= mz_shift_next_peptide)
          {
            mz_shifts_actual[mz_position] = std::numeric_limits<double>::quiet_NaN();
            mz_shifts_actual_indices[mz_position] = -1;
          }
        }
      }
    }

    // remove blacklisted peaks
    // loop over isotopes in peptides
    for (int isotope = 0; isotope < peaks_per_peptide_max_; ++isotope)
    {
      // loop over peptides
      for (int peptide = 0; peptide < (int) pattern.getMassShiftCount(); ++peptide)
      {
        int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1; // index in m/z shift list
        int peak_index = mz_shifts_actual_indices[mz_position]; // index of the peak in the spectrum
        if (peak_index != -1)
        {
          bool black = blacklist_[spectrum][peak_index].black;
          bool black_exception = blacklist_[spectrum][peak_index].black_exception_mass_shift_index == pattern.getMassShiftIndex() &&
                                 blacklist_[spectrum][peak_index].black_exception_charge == pattern.getCharge() &&
                                 blacklist_[spectrum][peak_index].black_exception_mz_position == mz_position;
          if (black && !black_exception)
          {
            mz_shifts_actual[mz_position] = std::numeric_limits<double>::quiet_NaN();
            mz_shifts_actual_indices[mz_position] = -1;
          }
        }
      }
    }

    // count how many isotopic peaks seen simultaneously in all of the peptides
    // and (optionally) remove peaks following missing ones
    int peaks_found_in_all_peptides = peaks_per_peptide_max_;
    for (int peptide = 0; peptide < (int) pattern.getMassShiftCount(); ++peptide)
    {
      bool missing_peak_seen = false;
      for (int isotope = 0; isotope < peaks_per_peptide_max_; ++isotope)
      {
        int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1; // index in m/z shift list
        int index = mz_shifts_actual_indices[mz_position]; // peak index in spectrum
        if (index == -1)
        {
          missing_peak_seen = true;
          peaks_found_in_all_peptides = min(peaks_found_in_all_peptides, isotope);
        }
        // if missing peaks are not allowed and we have already encountered one,
        // delete all higher isotopic peaks
        if (missing_peaks_ && missing_peak_seen)
        {
          mz_shifts_actual[mz_position] = std::numeric_limits<double>::quiet_NaN();
          mz_shifts_actual_indices[mz_position] = -1;
        }
      }
    }

    return peaks_found_in_all_peptides;
  }

  bool MultiplexFiltering::monoIsotopicPeakIntensityFilter(MultiplexPeakPattern pattern, int spectrum_index, const vector<int>& mz_shifts_actual_indices) const
  {
    MSExperiment<Peak1D>::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      int peak_index = mz_shifts_actual_indices[peptide * (peaks_per_peptide_max_ + 1) + 1];
      if (peak_index == -1)
      {
        // peak not found
        return true;
      }
      MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin() + peak_index;
      if (it_mz->getIntensity() < intensity_cutoff_)
      {
        // below intensity threshold
        return true;
      }
    }
    return false;
  }

  bool MultiplexFiltering::zerothPeakFilter(MultiplexPeakPattern pattern, const vector<double>& intensities_actual) const
  {
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // scaling factor for the zeroth peak intensity
      // (The zeroth peak is problematic if its intensity exceeds zero_scaling * intensity of mono-isotopic peak.)
      double zero_scaling = 0.7;
      if (boost::math::isnan(intensities_actual[peptide * (peaks_per_peptide_max_ + 1)]))
      {
        // zeroth peak not found
        continue;
      }
      else if (boost::math::isnan(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + 1]))
      {
        // first peak not found
        return true;
      }
      else if (intensities_actual[peptide * (peaks_per_peptide_max_ + 1)] > zero_scaling * intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + 1])
      {
        return true;
      }
    }

    return false;
  }

  bool MultiplexFiltering::peptideSimilarityFilter(MultiplexPeakPattern pattern, const vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline) const
  {
    std::vector<double> isotope_pattern_1;
    std::vector<double> isotope_pattern_2;
    for (unsigned peptide = 1; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      for (int isotope = 0; isotope < peaks_found_in_all_peptides_spline; ++isotope)
      {
        if (boost::math::isnan(intensities_actual[isotope + 1]))
        {
          // no peak found, hence assume the intensity at this position to be zero
          isotope_pattern_1.push_back(0);
        }
        else
        {
          isotope_pattern_1.push_back(intensities_actual[isotope + 1]);
        }
        if (boost::math::isnan(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1]))
        {
          // no peak found, hence assume the intensity at this position to be zero
          isotope_pattern_2.push_back(0);
        }
        else
        {
          isotope_pattern_2.push_back(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1]);
        }
      }
      if (getPatternSimilarity(isotope_pattern_1, isotope_pattern_2) < peptide_similarity_)
      {
        return false;
      }
    }

    return true;
  }

  bool MultiplexFiltering::averagineSimilarityFilter(MultiplexPeakPattern pattern, const vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline, double mz) const
  {
    // Use a more restrictive averagine similarity when we are searching for peptide singlets.
    double similarity;
    if (pattern.getMassShiftCount() == 1)
    {
      // We are detecting peptide singlets.
      similarity = averagine_similarity_ + averagine_similarity_scaling_*(1 - averagine_similarity_);
    }
    else
    {
      // We are detecting peptide doublets or triplets or ...
      similarity = averagine_similarity_;
    }
    
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      vector<double> isotope_pattern;
      for (int isotope = 0; isotope < peaks_found_in_all_peptides_spline; ++isotope)
      {
        if (boost::math::isnan(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1]))
        {
          // no peak found, hence assume the intensity at this position to be zero
          isotope_pattern.push_back(0);
        }
        else
        {
          isotope_pattern.push_back(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1]);
        }
      }
      if (getAveragineSimilarity(isotope_pattern, mz * pattern.getCharge()) < similarity)
      {
        return false;
      }
    }

    return true;
  }

  void MultiplexFiltering::blacklistPeaks(MultiplexPeakPattern pattern, int spectrum, const vector<int>& mz_shifts_actual_indices, int peaks_found_in_all_peptides_spline)
  {
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      for (int isotope = 0; isotope < peaks_found_in_all_peptides_spline; ++isotope)
      {
        int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1; // index in m/z shift list
        int peak_index;

        // blacklist peaks in this spectrum
        peak_index = mz_shifts_actual_indices[mz_position];
        if (peak_index != -1 && !blacklist_[spectrum][peak_index].black)
        {
          blacklist_[spectrum][peak_index].black = true;
          blacklist_[spectrum][peak_index].black_exception_mass_shift_index = pattern.getMassShiftIndex();
          blacklist_[spectrum][peak_index].black_exception_charge = pattern.getCharge();
          blacklist_[spectrum][peak_index].black_exception_mz_position = mz_position;
        }

        // blacklist peaks in previous spectrum
        peak_index = registry_[spectrum][mz_shifts_actual_indices[mz_position]].index_in_previous_spectrum;
        if (peak_index != -1 && !blacklist_[spectrum - 1][peak_index].black)
        {
          blacklist_[spectrum - 1][peak_index].black = true;
          blacklist_[spectrum - 1][peak_index].black_exception_mass_shift_index = pattern.getMassShiftIndex();
          blacklist_[spectrum - 1][peak_index].black_exception_charge = pattern.getCharge();
          blacklist_[spectrum - 1][peak_index].black_exception_mz_position = mz_position;
        }
        
        // blacklist peaks in spectrum before previous one
        if (peak_index != -1 && spectrum > 1)
        {
          int peak_index_2 = registry_[spectrum - 1][peak_index].index_in_previous_spectrum;
          if (peak_index_2 != -1 && !blacklist_[spectrum - 2][peak_index_2].black)
          {
            blacklist_[spectrum - 2][peak_index_2].black = true;
            blacklist_[spectrum - 2][peak_index_2].black_exception_mass_shift_index = pattern.getMassShiftIndex();
            blacklist_[spectrum - 2][peak_index_2].black_exception_charge = pattern.getCharge();
            blacklist_[spectrum - 2][peak_index_2].black_exception_mz_position = mz_position;
          }
        }
        
        // blacklist peaks in next spectrum
        peak_index = registry_[spectrum][mz_shifts_actual_indices[mz_position]].index_in_next_spectrum;
        if (peak_index != -1 && !blacklist_[spectrum + 1][peak_index].black)
        {
          blacklist_[spectrum + 1][peak_index].black = true;
          blacklist_[spectrum + 1][peak_index].black_exception_mass_shift_index = pattern.getMassShiftIndex();
          blacklist_[spectrum + 1][peak_index].black_exception_charge = pattern.getCharge();
          blacklist_[spectrum + 1][peak_index].black_exception_mz_position = mz_position;
        }
        
        // blacklist peaks in spectrum after next one
        if (peak_index != -1 && spectrum + 2 < (int) blacklist_.size())
        {
          int peak_index_2 = registry_[spectrum + 1][peak_index].index_in_next_spectrum;
          if (peak_index_2 != -1 && !blacklist_[spectrum + 2][peak_index_2].black)
          {
            blacklist_[spectrum + 2][peak_index_2].black = true;
            blacklist_[spectrum + 2][peak_index_2].black_exception_mass_shift_index = pattern.getMassShiftIndex();
            blacklist_[spectrum + 2][peak_index_2].black_exception_charge = pattern.getCharge();
            blacklist_[spectrum + 2][peak_index_2].black_exception_mz_position = mz_position;
          }
        }

      }
    }
  }

  int MultiplexFiltering::getPeakIndex(std::vector<double> peak_position, int start, double mz, double scaling) const
  {
    vector<int> valid_index; // indices of valid peaks that lie within the ppm range of the expected peak
    vector<double> valid_deviation; // ppm deviations between expected and (valid) actual peaks

    if (peak_position[start] < mz)
    {
      for (unsigned i = start; i < peak_position.size(); ++i)
      {
        double mz_min;
        double mz_max;
        if (mz_tolerance_unit_)
        {
          mz_min = (1 - scaling * mz_tolerance_ / 1000000) * peak_position[i];
          mz_max = (1 + scaling * mz_tolerance_ / 1000000) * peak_position[i];
        }
        else
        {
          mz_min = peak_position[i] - scaling * mz_tolerance_;
          mz_max = peak_position[i] + scaling * mz_tolerance_;
        }

        if (mz >= mz_min && mz <= mz_max)
        {
          valid_index.push_back(i);
          valid_deviation.push_back(abs(mz - peak_position[i]) / mz * 1000000);
        }
        if (mz < peak_position[i])
        {
          break;
        }
      }
    }
    else
    {
      for (int i = start; i >= 0; --i)
      {
        double mz_min;
        double mz_max;
        if (mz_tolerance_unit_)
        {
          mz_min = (1 - scaling * mz_tolerance_ / 1000000) * peak_position[i];
          mz_max = (1 + scaling * mz_tolerance_ / 1000000) * peak_position[i];
        }
        else
        {
          mz_min = peak_position[i] - scaling * mz_tolerance_;
          mz_max = peak_position[i] + scaling * mz_tolerance_;
        }

        if (mz >= mz_min && mz <= mz_max)
        {
          valid_index.push_back(i);
          valid_deviation.push_back(abs(mz - peak_position[i]) / mz * 1000000);
        }
        if (mz > peak_position[i])
        {
          break;
        }
      }
    }

    if (valid_index.size() == 0)
    {
      return -1;
    }
    else
    {
      // find best index
      int best_index = valid_index[0];
      double best_deviation = valid_deviation[0];
      for (unsigned i = 1; i < valid_index.size(); ++i)
      {
        if (valid_deviation[i] < best_deviation)
        {
          best_index = valid_index[i];
          best_deviation = valid_deviation[i];
        }
      }

      return best_index;
    }
  }

  double MultiplexFiltering::getPatternSimilarity(vector<double> pattern1, vector<double> pattern2) const
  {
    if (pattern1.empty() || pattern2.empty())
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    return OpenMS::Math::pearsonCorrelationCoefficient(pattern1.begin(), pattern1.end(), pattern2.begin(), pattern2.end());
  }

  double MultiplexFiltering::getAveragineSimilarity(vector<double> pattern, double m) const
  {
    // construct averagine distribution
    IsotopeDistribution distribution;
    vector<double> averagine_pattern;
    distribution.setMaxIsotope(pattern.size());
    distribution.estimateFromPeptideWeight(m);
    for (IsotopeDistribution::Iterator it = distribution.begin(); it != distribution.end(); ++it)
    {
      averagine_pattern.push_back(it->second);
    }

    return getPatternSimilarity(pattern, averagine_pattern);
  }

}
