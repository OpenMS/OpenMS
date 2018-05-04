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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <QDir>

using namespace std;
using namespace boost::math;

namespace OpenMS
{

  MultiplexFiltering::MultiplexFiltering(const PeakMap& exp_picked, const std::vector<MultiplexIsotopicPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averigine_type) :
    exp_picked_(exp_picked), patterns_(patterns), peaks_per_peptide_min_(peaks_per_peptide_min), peaks_per_peptide_max_(peaks_per_peptide_max), missing_peaks_(missing_peaks), intensity_cutoff_(intensity_cutoff), mz_tolerance_(mz_tolerance), mz_tolerance_unit_(mz_tolerance_unit), peptide_similarity_(peptide_similarity), averagine_similarity_(averagine_similarity), averagine_similarity_scaling_(averagine_similarity_scaling), averagine_type_(averigine_type)
  {
  }

  int MultiplexFiltering::positionsAndBlacklistFilter_(const MultiplexIsotopicPeakPattern& pattern, int spectrum,
                                                      const vector<double>& peak_position, int peak,
                                                      vector<double>& mz_shifts_actual,
                                                      vector<int>& mz_shifts_actual_indices) const
  {
    // Try to find peaks at the expected m/z positions
    // loop over expected m/z shifts of a peak pattern
    unsigned found_peaks(0);
    for (unsigned mz_position = 0; mz_position < pattern.getMZShiftCount(); ++mz_position)
    {
      double scaling = 1;
      if (mz_position % (peaks_per_peptide_max_ + 1) == 0)
      {
        // Let us be more lenient when looking for zeroths peaks
        // i.e. allow for an increased deviation between expected m/z position and the actual one
        scaling = 2;
      }

      int index = getPeakIndex_(peak_position, peak, peak_position[peak] + pattern.getMZShiftAt(mz_position), scaling);
      if (index != -1)
      {
        ++found_peaks;
        mz_shifts_actual.push_back(peak_position[index] - peak_position[peak]);
        mz_shifts_actual_indices.push_back(index);
      }
      else
      {
        mz_shifts_actual.push_back(std::numeric_limits<double>::quiet_NaN());
        mz_shifts_actual_indices.push_back(-1);
      }
    }

    // early out: Need to find at least (peaks_per_peptide * number_of_peptides) isotopic peaks.
    if (found_peaks < peaks_per_peptide_min_ * pattern.getMassShiftCount()) return -1;

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

  bool MultiplexFiltering::monoIsotopicPeakIntensityFilter_(const MultiplexIsotopicPeakPattern& pattern, int spectrum_index, const vector<int>& mz_shifts_actual_indices) const
  {
    PeakMap::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      int peak_index = mz_shifts_actual_indices[peptide * (peaks_per_peptide_max_ + 1) + 1];
      if (peak_index == -1)
      {
        // peak not found
        return true;
      }
      MSSpectrum::ConstIterator it_mz = it_rt->begin() + peak_index;
      if (it_mz->getIntensity() < intensity_cutoff_)
      {
        // below intensity threshold
        return true;
      }
    }
    return false;
  }

  bool MultiplexFiltering::zerothPeakFilter_(const MultiplexIsotopicPeakPattern& pattern, const vector<double>& intensities_actual) const
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

  bool MultiplexFiltering::peptideSimilarityFilter_(const MultiplexIsotopicPeakPattern& pattern, const vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline) const
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
      if (getPatternSimilarity_(isotope_pattern_1, isotope_pattern_2) < peptide_similarity_)
      {
        return false;
      }
    }

    return true;
  }

  bool MultiplexFiltering::averagineSimilarityFilter_(const MultiplexIsotopicPeakPattern& pattern, const vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline, double mz) const
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
      if (getAveragineSimilarity_(isotope_pattern, mz * pattern.getCharge()) < similarity)
      {
        return false;
      }
    }

    return true;
  }

  void MultiplexFiltering::blacklistPeaks_(const MultiplexIsotopicPeakPattern& pattern, int spectrum, const vector<int>& mz_shifts_actual_indices, int peaks_found_in_all_peptides_spline)
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

  int MultiplexFiltering::getPeakIndex_(const std::vector<double>& peak_position, int start, double mz, double scaling) const
  {
    const double tolerance_th = mz_tolerance_unit_ ? (scaling * mz_tolerance_ / 1000000) * peak_position[start] : scaling * mz_tolerance_;
    const double mz_min = mz - tolerance_th;
    const double mz_max = mz + tolerance_th;

    std::vector<double>::const_iterator lb = std::lower_bound(peak_position.begin(), peak_position.end(), mz_min);
    std::vector<double>::const_iterator ub = std::upper_bound(lb, peak_position.end(), mz_max);

    double smallest_error = scaling * mz_tolerance_; // initialize to the maximum  allowed error 
    int smallest_error_index = -1;

    for (; lb != ub; ++lb)
    {
      const double error = abs(*lb - mz);
      if (error <= smallest_error)
      {
        smallest_error = error;
        smallest_error_index = lb - peak_position.begin();
      }
    }    

    return smallest_error_index;
  }

  double MultiplexFiltering::getPatternSimilarity_(const vector<double>& pattern1, const vector<double>& pattern2) const
  {
    if (pattern1.empty() || pattern2.empty())
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    return OpenMS::Math::pearsonCorrelationCoefficient(pattern1.begin(), pattern1.end(), pattern2.begin(), pattern2.end());
  }


  double MultiplexFiltering::getAveragineSimilarity_(const vector<double>& pattern, double m) const

  {
    // construct averagine distribution
    CoarseIsotopePatternGenerator solver(pattern.size());
    IsotopeDistribution distribution;
    vector<double> averagine_pattern;
    if (averagine_type_ == "peptide")
    {
        distribution = solver.estimateFromPeptideWeight(m);
    }
    else if (averagine_type_ == "RNA")
    {
      distribution = solver.estimateFromRNAWeight(m);
    }
    else if (averagine_type_ == "DNA")
    {
        distribution = solver.estimateFromDNAWeight(m);
    }
    else
    {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Averagine type unrecognized.");;
    }

    for (IsotopeDistribution::Iterator it = distribution.begin(); it != distribution.end(); ++it)
    {
      averagine_pattern.push_back(it->getIntensity());
    }

    return getPatternSimilarity_(pattern, averagine_pattern);
  }

}
