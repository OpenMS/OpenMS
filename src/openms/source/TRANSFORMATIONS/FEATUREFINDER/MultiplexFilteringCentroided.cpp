// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
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

  MultiplexFilteringCentroided::MultiplexFilteringCentroided(MSExperiment<Peak1D> exp_picked, std::vector<MultiplexPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String out_debug) :
    MultiplexFiltering(exp_picked, patterns, peaks_per_peptide_min, peaks_per_peptide_max, missing_peaks, intensity_cutoff, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, out_debug)
  {

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
          reference.index_in_previous_spectrum = getPeakIndex(index - 1, it_mz->getMZ(), 3.0);
        }
        else
        {
          reference.index_in_previous_spectrum = -1;
        }
        if (index + 1 < (int) exp_picked_.size())
        {
          reference.index_in_next_spectrum = getPeakIndex(index + 1, it_mz->getMZ(), 3.0);
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

  vector<MultiplexFilterResult> MultiplexFilteringCentroided::filter()
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns_.size() * exp_picked_.size(), "filtering LC-MS data");

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
      for (MSExperiment<Peak1D>::Iterator it_rt_picked = exp_picked_.begin(); it_rt_picked < exp_picked_.end(); ++it_rt_picked)
      {
        setProgress(++progress);

        int spectrum = it_rt_picked - exp_picked_.begin(); // index of the spectrum in exp_picked_
        double rt_picked = it_rt_picked->getRT();

        // vectors of peak details
        vector<double> peak_position;
        vector<double> peak_intensity;
        for (MSSpectrum<Peak1D>::Iterator it_mz = it_rt_picked->begin(); it_mz < it_rt_picked->end(); ++it_mz)
        {
          peak_position.push_back(it_mz->getMZ());
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
            if (debug_)
            {
              Peak2D data_point;
              data_point.setRT(rt_picked);
              data_point.setMZ(peak_position[peak]);
              data_point.setIntensity(1); // filter 1 failed
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
              data_point.setIntensity(2); // filter 2 failed
              debug_rejected.push_back(data_point);
            }
            continue;
          }

          /**
           * Filter (3): non-local intensity filter
           * Are the peak intensities of all peptides above the cutoff?
           */
          std::vector<double> intensities_actual; // peak intensities @ m/z peak position + actual m/z shift
          int peaks_found_in_all_peptides_centroided = nonLocalIntensityFilter(patterns_[pattern], spectrum, mz_shifts_actual_indices, intensities_actual, peaks_found_in_all_peptides);
          if (peaks_found_in_all_peptides_centroided < peaks_per_peptide_min_)
          {
            if (debug_)
            {
              Peak2D data_point;
              data_point.setRT(rt_picked);
              data_point.setMZ(peak_position[peak]);
              data_point.setIntensity(3); // filter 3 failed
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
              data_point.setMZ(peak_position[peak]);
              data_point.setIntensity(4); // filter 4 failed
              debug_rejected.push_back(data_point);
            }
            continue;
          }

          /**
           * Filter (5): peptide similarity filter
           * How similar are the isotope patterns of the peptides?
           */
          bool peptide_similarity = peptideSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_centroided);
          if (!peptide_similarity)
          {
            if (debug_)
            {
              Peak2D data_point;
              data_point.setRT(rt_picked);
              data_point.setMZ(peak_position[peak]);
              data_point.setIntensity(5); // filter 5 failed
              debug_rejected.push_back(data_point);
            }
            continue;
          }

          /**
           * Filter (6): averagine similarity filter
           * Does each individual isotope pattern resemble a peptide?
           */
          bool averagine_similarity = averagineSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_centroided, peak_position[peak]);
          if (!averagine_similarity)
          {
            if (debug_)
            {
              Peak2D data_point;
              data_point.setRT(rt_picked);
              data_point.setMZ(peak_position[peak]);
              data_point.setIntensity(6); // filter 6 failed
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
            data_point.setMZ(peak_position[peak]);
            data_point.setIntensity(intensities_actual[1]); // all filters passed
            debug_filtered.push_back(data_point);
          }

          // add the peak to the result
          vector<MultiplexFilterResultRaw> results_raw;
          result.addFilterResultPeak(peak_position[peak], rt_picked, mz_shifts_actual, intensities_actual, results_raw);

          // blacklist peaks in the current spectrum and the two neighbouring ones
          blacklistPeaks(patterns_[pattern], spectrum, mz_shifts_actual_indices, peaks_found_in_all_peptides_centroided);

        }
      }

      // add results of this pattern to list
      filter_results.push_back(result);

      // write debug output
      if (debug_)
      {
        // Writes for each peak pattern two debug files.
        // One containing rejected data points. The intensity encodes which of the six filters failed.
        // The second file containing data points that passed all filters.
        writeDebug(pattern, true, debug_rejected);
        writeDebug(pattern, false, debug_filtered);
      }
    }

    endProgress();

    return filter_results;
  }

  int MultiplexFilteringCentroided::nonLocalIntensityFilter(MultiplexPeakPattern pattern, int spectrum_index, const std::vector<int>& mz_shifts_actual_indices, std::vector<double>& intensities_actual, int peaks_found_in_all_peptides) const
  {
    MSExperiment<Peak1D>::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;

    // read out peak intensities
    for (int i = 0; i < (int) mz_shifts_actual_indices.size(); ++i)
    {
      int peak_index = mz_shifts_actual_indices[i];
      if (peak_index != -1)
      {
        intensities_actual.push_back((it_rt->begin() + peak_index)->getIntensity());
      }
      else
      {
        intensities_actual.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }

    for (int isotope = 0; isotope < peaks_found_in_all_peptides; ++isotope)
    {
      bool seen_in_all_peptides = true;
      for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
      {
        int peak_index = mz_shifts_actual_indices[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1];
        if (peak_index == -1)
        {
          // peak not found
          seen_in_all_peptides = false;
          break;
        }
        else if ((it_rt->begin() + peak_index)->getIntensity() < intensity_cutoff_)
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

  int MultiplexFilteringCentroided::getPeakIndex(int spectrum_index, double mz, double scaling) const
  {
    MSExperiment<Peak1D>::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;
    MSSpectrum<Peak1D>::ConstIterator it_mz;
    for (it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
    {
      double mz_min;
      double mz_max;
      if (mz_tolerance_unit_)
      {
        mz_min = (1 - scaling * mz_tolerance_ / 1000000) * it_mz->getMZ();
        mz_max = (1 + scaling * mz_tolerance_ / 1000000) * it_mz->getMZ();
      }
      else
      {
        mz_min = it_mz->getMZ() - scaling * mz_tolerance_;
        mz_max = it_mz->getMZ() + scaling * mz_tolerance_;
      }

      if (mz >= mz_min && mz <= mz_max)
      {
        return it_mz - it_rt->begin();
      }

      if (mz < mz_min)
      {
        return -1;
      }
    }

    return -1;
  }

}
