// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include<QDir>

using namespace std;
using namespace boost::math;

namespace OpenMS
{

  MultiplexFiltering::MultiplexFiltering(MSExperiment<Peak1D> exp_profile, MSExperiment<Peak1D> exp_picked, vector<vector<PeakPickerHiRes::PeakBoundary> > boundaries, std::vector<MultiplexPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, String out_debug) :
    exp_profile_(exp_profile), exp_picked_(exp_picked), boundaries_(boundaries), patterns_(patterns), peaks_per_peptide_min_(peaks_per_peptide_min), peaks_per_peptide_max_(peaks_per_peptide_max), missing_peaks_(missing_peaks), intensity_cutoff_(intensity_cutoff), mz_tolerance_(mz_tolerance), mz_tolerance_unit_(mz_tolerance_unit), peptide_similarity_(peptide_similarity), averagine_similarity_(averagine_similarity), out_debug_(out_debug), debug_(out_debug.trim().length() > 0)
  {
    if (exp_profile_.size() != exp_picked_.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Profile and centroided data do not contain same number of spectra.");
    }

    if (exp_picked_.size() != boundaries_.size())
    {
      stringstream stream;
      stream << "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra. (";
      stream << exp_picked_.size();
      stream << "!=";
      stream << boundaries_.size();
      stream << ")";
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream.str());
    }

    // fill peak registry and initialise blacklist
    MSExperiment<Peak1D>::Iterator it_rt;
    vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
    for (it_rt = exp_picked_.begin(), it_rt_boundaries = boundaries_.begin();
         it_rt < exp_picked_.end() && it_rt_boundaries < boundaries_.end();
         ++it_rt, ++it_rt_boundaries)
    {
      int index = it_rt - exp_picked_.begin();

      vector<PeakReference> registry_spec;
      vector<BlackListEntry> blacklist_spec;
      for (MSSpectrum<Peak1D>::Iterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
      {
        // peak registry
        PeakReference reference;
        if (index > 0)
        {
          reference.index_in_previous_spectrum = getPeakIndex(index - 1, it_mz->getMZ(), 1.0);
        }
        else
        {
          reference.index_in_previous_spectrum = -1;
        }
        if (index + 1 < (int) exp_picked_.size())
        {
          reference.index_in_next_spectrum = getPeakIndex(index + 1, it_mz->getMZ(), 1.0);
        }
        else
        {
          reference.index_in_next_spectrum = -1;
        }
        registry_spec.push_back(reference);

        // blacklist
        BlackListEntry entry;
        entry.black = false;
        entry.black_exception_mass_shift_index = -1;
        entry.black_exception_charge = -1;
        entry.black_exception_mz_position = -1;
        blacklist_spec.push_back(entry);
      }
      registry_.push_back(registry_spec);
      blacklist_.push_back(blacklist_spec);
    }

  }

  vector<MultiplexFilterResult> MultiplexFiltering::filter()
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns_.size() * exp_profile_.size(), "filtering LC-MS data");
      
    // list of filter results for each peak pattern
    vector<MultiplexFilterResult> filter_results;

    // loop over patterns
    for (unsigned pattern = 0; pattern < patterns_.size(); ++pattern)
    {
      // data structure storing peaks which pass all filters
      MultiplexFilterResult result;

      // m/z position is rejected by a particular filter (or passing all of them)
      vector<Peak2D> debug_rejected;
      vector<Peak2D> debug_filtered;

      // iterate over spectra
      MSExperiment<Peak1D>::Iterator it_rt_profile;
      MSExperiment<Peak1D>::Iterator it_rt_picked;
      vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
      for (it_rt_profile = exp_profile_.begin(), it_rt_picked = exp_picked_.begin(), it_rt_boundaries = boundaries_.begin();
           it_rt_profile < exp_profile_.end() && it_rt_picked < exp_picked_.end() && it_rt_boundaries < boundaries_.end();
           ++it_rt_profile, ++it_rt_picked, ++it_rt_boundaries)
      {
        if ((*it_rt_picked).size() != (*it_rt_boundaries).size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Number of peaks and number of peak boundaries differ.");
        }
        
        setProgress(++progress);

        int spectrum = it_rt_profile - exp_profile_.begin();            // index of the spectrum in exp_profile_, exp_picked_ and boundaries_
        double rt_picked = it_rt_picked->getRT();
        
        // spline fit profile data
        SplineSpectrum spline(*it_rt_profile);
        SplineSpectrum::Navigator nav = spline.getNavigator();

        // vectors of peak details
        vector<double> peak_position;
        vector<double> peak_min;
        vector<double> peak_max;
        vector<double> peak_intensity;
        MSSpectrum<Peak1D>::Iterator it_mz;
        vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary;
        for (it_mz = it_rt_picked->begin(), it_mz_boundary = it_rt_boundaries->begin();
             it_mz < it_rt_picked->end(), it_mz_boundary < it_rt_boundaries->end();
             ++it_mz, ++it_mz_boundary)
        {
          peak_position.push_back(it_mz->getMZ());
          peak_min.push_back((*it_mz_boundary).mz_min);
          peak_max.push_back((*it_mz_boundary).mz_max);
          peak_intensity.push_back(it_mz->getIntensity());
        }

        // iterate over peaks in spectrum (mz)
        for (unsigned peak = 0; peak < peak_position.size(); ++peak)
        {

          /**
           * Filter (1): m/z position and blacklist filter
           * Are there non-black peaks with the expected relative m/z shifts?
           */
          vector<double> mz_shifts_actual;              // actual m/z shifts (differ slightly from expected m/z shifts)
          vector<int> mz_shifts_actual_indices;              // peak indices in the spectrum corresponding to the actual m/z shifts
          int peaks_found_in_all_peptides = positionsAndBlacklistFilter(patterns_[pattern], spectrum, peak_position, peak, mz_shifts_actual, mz_shifts_actual_indices);
          if (peaks_found_in_all_peptides < peaks_per_peptide_min_)
          {
            if (debug_)
            {
              Peak2D data_point;
              data_point.setRT(rt_picked);
              data_point.setMZ(peak_position[peak]);
              data_point.setIntensity(1);                  // filter 1 failed
              debug_rejected.push_back(data_point);
            }
            continue;
          }

          /**
            * Filter (2): blunt intensity filter
            * Are the mono-isotopic peak intensities of all peptides above the cutoff?
            */
          bool bluntVeto = monoIsotopicPeakIntensityFilter(patterns_[pattern], spectrum, mz_shifts_actual_indices);
          if (bluntVeto)
          {
            if (debug_)
            {
              Peak2D data_point;
              data_point.setRT(rt_picked);
              data_point.setMZ(peak_position[peak]);
              data_point.setIntensity(2);                  // filter 2 failed
              debug_rejected.push_back(data_point);
            }
            continue;
          }

          // Arrangement of peaks looks promising. Now scan through the spline fitted data.
          vector<MultiplexFilterResultRaw> results_raw;              // raw data points of this peak that will pass the remaining filters
          bool blacklisted = false;              // Has this peak already been blacklisted?
          for (double mz = peak_min[peak]; mz < peak_max[peak]; mz = nav.getNextMz(mz))
          {
            /**
             * Filter (3): non-local intensity filter
             * Are the spline interpolated intensities at m/z above the threshold?
             */
            vector<double> intensities_actual;                // spline interpolated intensities @ m/z + actual m/z shift
            int peaks_found_in_all_peptides_spline = nonLocalIntensityFilter(patterns_[pattern], mz_shifts_actual, mz_shifts_actual_indices, nav, intensities_actual, peaks_found_in_all_peptides, mz);
            if (peaks_found_in_all_peptides_spline < peaks_per_peptide_min_)
            {
              if (debug_)
              {
                Peak2D data_point;
                data_point.setRT(rt_picked);
                data_point.setMZ(mz);
                data_point.setIntensity(3);                    // filter 3 failed
                debug_rejected.push_back(data_point);
              }
              continue;
            }

            /**
             * Filter (4): zeroth peak filter
             * There should not be a significant peak to the left of the mono-isotopic
             * (i.e. first) peak.
             */
            bool zero_peak = zerothPeakFilter(patterns_[pattern], intensities_actual);
            if (zero_peak)
            {
              if (debug_)
              {
                Peak2D data_point;
                data_point.setRT(rt_picked);
                data_point.setMZ(mz);
                data_point.setIntensity(4);                    // filter 4 failed
                debug_rejected.push_back(data_point);
              }
              continue;
            }

            /**
             * Filter (5): peptide similarity filter
             * How similar are the isotope patterns of the peptides?
             */
            vector<double> isotope_pattern_1;
            vector<double> isotope_pattern_2;
            bool peptide_similarity = peptideSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_spline, isotope_pattern_1, isotope_pattern_2);
            if (!peptide_similarity)
            {
              if (debug_)
              {
                Peak2D data_point;
                data_point.setRT(rt_picked);
                data_point.setMZ(mz);
                data_point.setIntensity(5);                    // filter 5 failed
                debug_rejected.push_back(data_point);
              }
              continue;
            }

            /**
             * Filter (6): averagine similarity filter
             * Does each individual isotope pattern resemble a peptide?
             */
            bool averagine_similarity = averagineSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_spline, mz);
            if (!averagine_similarity)
            {
              if (debug_)
              {
                Peak2D data_point;
                data_point.setRT(rt_picked);
                data_point.setMZ(mz);
                data_point.setIntensity(6);                    // filter 6 failed
                debug_rejected.push_back(data_point);
              }
              continue;
            }

            /**
             * All filters passed.
             */
            if (debug_)
            {
              Peak2D data_point;
              data_point.setRT(rt_picked);
              data_point.setMZ(mz);
              data_point.setIntensity(intensities_actual[1]);                  // all filters passed
              debug_filtered.push_back(data_point);
            }

            // add raw data point to list that passed all filters
            MultiplexFilterResultRaw result_raw(mz, mz_shifts_actual, intensities_actual);
            results_raw.push_back(result_raw);

            // blacklist peaks in the current spectrum and the two neighbouring ones
            if (!blacklisted)
            {
              blacklistPeaks(patterns_[pattern], spectrum, mz_shifts_actual_indices, peaks_found_in_all_peptides_spline);
              blacklisted = true;
            }

          }

          // add the peak with its corresponding raw data to the result
          if (results_raw.size() > 2)
          {
            // Scanning over the profile of the peak, we want at least three raw data points to pass all filters.
            vector<double> intensities_actual;
            for (unsigned i = 0; i < mz_shifts_actual_indices.size(); ++i)
            {
              int index = mz_shifts_actual_indices[i];
              if (index == -1)
              {
                // no peak found
                intensities_actual.push_back(std::numeric_limits<double>::quiet_NaN());
              }
              else
              {
                intensities_actual.push_back(peak_intensity[mz_shifts_actual_indices[i]]);
              }
            }
            result.addFilterResultPeak(peak_position[peak], rt_picked, mz_shifts_actual, intensities_actual, results_raw);
          }

        }

      }

      // add results of this pattern to list
      filter_results.push_back(result);

      // write debug output
      if (debug_)
      {
        // Writes for each peak pattern two debug files.
        // One containg rejected data points. The intensity encodes which of the six filters failed.
        // The second file containing data points that passed all filters.
        writeDebug(pattern, true, debug_rejected);
        writeDebug(pattern, false, debug_filtered);
      }
    }

    endProgress();

    return filter_results;
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
      double mz_shift_next_peptide = mz_shifts_actual[(peptide + 1) * (peaks_per_peptide_max_ + 1) + 1];          // m/z shift of the mono-isotopic peak of the following peptide
      if (!(boost::math::isnan)(mz_shift_next_peptide))
      {
        for (int isotope = 0; isotope < peaks_per_peptide_max_; ++isotope)
        {
          int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1;              // index in m/z shift list
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
        int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1;            // index in m/z shift list
        int peak_index = mz_shifts_actual_indices[mz_position];            // index of the peak in the spectrum
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
        int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1;           // index in m/z shift list
        int index = mz_shifts_actual_indices[mz_position];            // peak index in spectrum
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

  int MultiplexFiltering::nonLocalIntensityFilter(MultiplexPeakPattern pattern, const vector<double>& mz_shifts_actual, const vector<int>& mz_shifts_actual_indices, SplineSpectrum::Navigator nav, std::vector<double>& intensities_actual, int peaks_found_in_all_peptides, double mz) const
  {
    // calculate intensities
    for (int i = 0; i < (int) mz_shifts_actual_indices.size(); ++i)
    {
      if (mz_shifts_actual_indices[i] != -1)
      {
        intensities_actual.push_back(nav.eval(mz + mz_shifts_actual[i]));
      }
      else
      {
        intensities_actual.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }

    // peaks_found_in_all_peptides is the number of peaks in this region. At this particular m/z
    // some of the intensities might be below the cutoff. => return isotope
    for (int isotope = 0; isotope < peaks_found_in_all_peptides; ++isotope)
    {
      bool seen_in_all_peptides = true;
      for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
      {
        if (boost::math::isnan(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1]))
        {
          // peak not found
          seen_in_all_peptides = false;
        }
        else if (intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1] < intensity_cutoff_)
        {
          // below intensity threshold
          seen_in_all_peptides = false;
        }
      }
      if (!seen_in_all_peptides)
      {
        return isotope;
      }
    }

    return peaks_found_in_all_peptides;
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

  bool MultiplexFiltering::peptideSimilarityFilter(MultiplexPeakPattern pattern, const vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline, vector<double> isotope_pattern_1, vector<double> isotope_pattern_2) const
  {
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
      if (getAveragineSimilarity(isotope_pattern, mz * pattern.getCharge()) < averagine_similarity_)
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
        int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1;            // index in m/z shift list
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

        // blacklist peaks in next spectrum
        peak_index = registry_[spectrum][mz_shifts_actual_indices[mz_position]].index_in_next_spectrum;
        if (peak_index != -1 && !blacklist_[spectrum + 1][peak_index].black)
        {
          blacklist_[spectrum + 1][peak_index].black = true;
          blacklist_[spectrum + 1][peak_index].black_exception_mass_shift_index = pattern.getMassShiftIndex();
          blacklist_[spectrum + 1][peak_index].black_exception_charge = pattern.getCharge();
          blacklist_[spectrum + 1][peak_index].black_exception_mz_position = mz_position;
        }

      }
    }
  }

  void MultiplexFiltering::writeDebug(int pattern, bool rejected, vector<Peak2D> points) const
  {
    MSExperiment<Peak1D> exp_debug;
    MSSpectrum<Peak1D> spec_debug;

    double rt = std::numeric_limits<double>::quiet_NaN();
    int spec_id = 0;
    for (vector<Peak2D>::const_iterator it = points.begin(); it != points.end(); ++it)
    {
      if ((boost::math::isnan)(rt) || (*it).getRT() > rt)
      {
        if (!(boost::math::isnan)(rt))
        {
          exp_debug.addSpectrum(spec_debug);
          ++spec_id;
        }

        rt = (*it).getRT();
        spec_debug.clear(true);
        spec_debug.setRT(rt);
        spec_debug.setMSLevel(1);
        spec_debug.setNativeID(String("spectrum = ") + spec_id);
      }

      Peak1D peak;
      peak.setMZ((*it).getMZ());
      peak.setIntensity((*it).getIntensity());
      spec_debug.push_back(peak);

    }

    MzMLFile fileSpline;
    String file_name;
    if (rejected)
    {
      file_name = "debug_rejected_";
    }
    else
    {
      file_name = "debug_passed_";
    }
    QDir dir(out_debug_.toQString());
    if (!dir.cdUp())
    {
        std::stringstream stream;
        stream << "Could not navigate to directory for debug output '" << String(dir.dirName()) << "'.";
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream.str());
    }
    if (!dir.exists() && !dir.mkpath("."))
    {
        std::stringstream stream;
        stream << "Could not create directory for debug output '" << String(dir.dirName()) << "'.";
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, stream.str());
    }
    file_name = out_debug_ + "/" + file_name + pattern + ".mzML";    // Correct way of writing to absolute path?
    fileSpline.store(file_name, exp_debug);

  }

  int MultiplexFiltering::getPeakIndex(int spectrum_index, double mz, double scaling) const
  {
    MSExperiment<Peak1D>::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;
    vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries = boundaries_.begin() + spectrum_index;

    MSSpectrum<Peak1D>::ConstIterator it_mz;
    vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundaries;
    for (it_mz = it_rt->begin(), it_mz_boundaries = it_rt_boundaries->begin();
         it_mz < it_rt->end(), it_mz_boundaries < it_rt_boundaries->end();
         ++it_mz, ++it_mz_boundaries)
    {
      if (mz >= scaling * (*it_mz_boundaries).mz_min + (1 - scaling) * it_mz->getMZ() &&
          mz <= scaling * (*it_mz_boundaries).mz_max + (1 - scaling) * it_mz->getMZ())
      {
        return it_mz - it_rt->begin();
      }
      if (mz < scaling * (*it_mz_boundaries).mz_min + (1 - scaling) * it_mz->getMZ())
      {
        return -1;
      }
    }

    return -1;
  }

  int MultiplexFiltering::getPeakIndex(std::vector<double> peak_position, int start, double mz, double scaling) const
  {
    vector<int> valid_index;        // indices of valid peaks that lie within the ppm range of the expected peak
    vector<double> valid_deviation;        // ppm deviations between expected and (valid) actual peaks

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
