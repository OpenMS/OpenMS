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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultPeak.h>
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

  MultiplexFilteringProfile::MultiplexFilteringProfile(const MSExperiment<Peak1D>& exp_profile, const MSExperiment<Peak1D>& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, const std::vector<MultiplexPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling) :
    MultiplexFiltering(exp_picked, patterns, peaks_per_peptide_min, peaks_per_peptide_max, missing_peaks, intensity_cutoff, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling), exp_profile_(exp_profile), boundaries_(boundaries)
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
    for (it_rt = exp_picked_.begin(); it_rt < exp_picked_.end(); ++it_rt)
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

  vector<MultiplexFilterResult> MultiplexFilteringProfile::filter()
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

      // iterate over spectra
      MSExperiment<Peak1D>::Iterator it_rt_profile;
      MSExperiment<Peak1D>::ConstIterator it_rt_picked;
      vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
      for (it_rt_profile = exp_profile_.begin(), it_rt_picked = exp_picked_.begin(), it_rt_boundaries = boundaries_.begin();
           it_rt_profile < exp_profile_.end() && it_rt_picked < exp_picked_.end() && it_rt_boundaries < boundaries_.end();
           ++it_rt_profile, ++it_rt_picked, ++it_rt_boundaries)
      {
        // skip empty spectra
        if ((*it_rt_profile).size() == 0 || (*it_rt_picked).size() == 0 || (*it_rt_boundaries).size() == 0)
        {
          continue;
        }
        
        if ((*it_rt_picked).size() != (*it_rt_boundaries).size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Number of peaks and number of peak boundaries differ.");
        }

        setProgress(++progress);

        int spectrum = it_rt_profile - exp_profile_.begin(); // index of the spectrum in exp_profile_, exp_picked_ and boundaries_
        double rt_picked = it_rt_picked->getRT();

        // spline fit profile data
        SplineSpectrum spline(*it_rt_profile);
        SplineSpectrum::Navigator nav = spline.getNavigator();

        // vectors of peak details
        vector<double> peak_position;
        vector<double> peak_min;
        vector<double> peak_max;
        vector<double> peak_intensity;
        MSSpectrum<Peak1D>::ConstIterator it_mz;
        vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary;
        for (it_mz = it_rt_picked->begin(), it_mz_boundary = it_rt_boundaries->begin();
             it_mz < it_rt_picked->end() && it_mz_boundary < it_rt_boundaries->end();
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
          vector<double> mz_shifts_actual; // actual m/z shifts (differ slightly from expected m/z shifts)
          vector<int> mz_shifts_actual_indices; // peak indices in the spectrum corresponding to the actual m/z shifts
          int peaks_found_in_all_peptides = positionsAndBlacklistFilter(patterns_[pattern], spectrum, peak_position, peak, mz_shifts_actual, mz_shifts_actual_indices);
          if (peaks_found_in_all_peptides < peaks_per_peptide_min_)
          {
            continue;
          }

          /**
           * Filter (2): blunt intensity filter
           * Are the mono-isotopic peak intensities of all peptides above the cutoff?
           */
          bool bluntVeto = monoIsotopicPeakIntensityFilter(patterns_[pattern], spectrum, mz_shifts_actual_indices);
          if (bluntVeto)
          {
            continue;
          }

          // Arrangement of peaks looks promising. Now scan through the spline fitted data.
          vector<MultiplexFilterResultRaw> results_raw; // raw data points of this peak that will pass the remaining filters
          bool blacklisted = false; // Has this peak already been blacklisted?
          for (double mz = peak_min[peak]; mz < peak_max[peak]; mz = nav.getNextMz(mz))
          {
            /**
             * Filter (3): non-local intensity filter
             * Are the spline interpolated intensities at m/z above the threshold?
             */
            vector<double> intensities_actual; // spline interpolated intensities @ m/z + actual m/z shift
            int peaks_found_in_all_peptides_spline = nonLocalIntensityFilter(patterns_[pattern], mz_shifts_actual, mz_shifts_actual_indices, nav, intensities_actual, peaks_found_in_all_peptides, mz);
            if (peaks_found_in_all_peptides_spline < peaks_per_peptide_min_)
            {
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
              continue;
            }

            /**
             * Filter (5): peptide similarity filter
             * How similar are the isotope patterns of the peptides?
             */
            bool peptide_similarity = peptideSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_spline);
            if (!peptide_similarity)
            {
              continue;
            }

            /**
             * Filter (6): averagine similarity filter
             * Does each individual isotope pattern resemble a peptide?
             */
            bool averagine_similarity = averagineSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_spline, mz);
            if (!averagine_similarity)
            {
              continue;
            }

            /**
             * All filters passed.
             */
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
    }

    endProgress();

    return filter_results;
  }

  int MultiplexFilteringProfile::nonLocalIntensityFilter(MultiplexPeakPattern pattern, const vector<double>& mz_shifts_actual, const vector<int>& mz_shifts_actual_indices, SplineSpectrum::Navigator nav, std::vector<double>& intensities_actual, int peaks_found_in_all_peptides, double mz) const
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
          break;
        }
        else if (intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1] < intensity_cutoff_)
        {
          // below intensity threshold
          seen_in_all_peptides = false;
          break;
        }
      }
      if (!seen_in_all_peptides)
      {
        return isotope;
      }
    }

    return peaks_found_in_all_peptides;
  }

  int MultiplexFilteringProfile::getPeakIndex(int spectrum_index, double mz, double scaling) const
  {
    MSExperiment<Peak1D>::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;
    vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries = boundaries_.begin() + spectrum_index;

    MSSpectrum<Peak1D>::ConstIterator it_mz;
    vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundaries;
    for (it_mz = it_rt->begin(), it_mz_boundaries = it_rt_boundaries->begin();
         it_mz < it_rt->end() && it_mz_boundaries < it_rt_boundaries->end();
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

}
