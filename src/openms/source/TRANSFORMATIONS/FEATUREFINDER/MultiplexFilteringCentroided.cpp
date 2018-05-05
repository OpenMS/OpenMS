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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <QDir>

using namespace std;
using namespace boost::math;

namespace OpenMS
{

  MultiplexFilteringCentroided::MultiplexFilteringCentroided(const PeakMap& exp_picked, const std::vector<MultiplexIsotopicPeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type) :
    MultiplexFiltering(exp_picked, patterns, peaks_per_peptide_min, peaks_per_peptide_max, missing_peaks, intensity_cutoff, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type)
  {

    // fill peak registry and initialise blacklist
    PeakMap::Iterator it_rt;
    blacklist_.reserve(exp_picked_.getNrSpectra());
    registry_.reserve(exp_picked_.getNrSpectra());
    for (it_rt = exp_picked_.begin(); it_rt < exp_picked_.end(); ++it_rt)
    {
      int index = it_rt - exp_picked_.begin();

      vector<PeakReference> registry_spec;
      vector<BlackListEntry> blacklist_spec;
      registry_spec.reserve(it_rt->size());
      blacklist_spec.reserve(it_rt->size());
      for (MSSpectrum::Iterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
      {
        const double scaling = 3.0; // mz_tolerance_ is the tolerance when searching for peak patterns. Let's be more generous when searching for corresponding peaks in neighbouring spectra.
        const double tolerance_th = mz_tolerance_unit_ ? (scaling * mz_tolerance_ / 1000000) * it_mz->getMZ() : scaling * mz_tolerance_;
        // peak registry
        PeakReference reference;
        if (index > 0)
        {
          reference.index_in_previous_spectrum = exp_picked_[index - 1].findNearest(it_mz->getMZ(), tolerance_th);
        }
        else
        {
          reference.index_in_previous_spectrum = -1;
        }
        if (index + 1 < (int) exp_picked_.size())
        {
          reference.index_in_next_spectrum = exp_picked_[index + 1].findNearest(it_mz->getMZ(), tolerance_th);
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

      // iterate over spectra
      for (PeakMap::Iterator it_rt_picked = exp_picked_.begin(); it_rt_picked < exp_picked_.end(); ++it_rt_picked)
      {
        // skip empty spectra
        if (it_rt_picked->empty())
        {
          continue;
        }

        setProgress(++progress);

        int spectrum = it_rt_picked - exp_picked_.begin(); // index of the spectrum in exp_picked_
        double rt_picked = it_rt_picked->getRT();

        // vectors of peak details
        vector<double> peak_position;
        vector<double> peak_intensity;
        peak_position.reserve(it_rt_picked->size());
        peak_intensity.reserve(it_rt_picked->size());
        for (MSSpectrum::Iterator it_mz = it_rt_picked->begin(); it_mz < it_rt_picked->end(); ++it_mz)
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

          mz_shifts_actual.reserve(patterns_[pattern].getMZShiftCount());
          mz_shifts_actual_indices.reserve(patterns_[pattern].getMZShiftCount());

          int peaks_found_in_all_peptides = positionsAndBlacklistFilter_(patterns_[pattern], spectrum, peak_position, peak, mz_shifts_actual, mz_shifts_actual_indices);
          if (peaks_found_in_all_peptides < peaks_per_peptide_min_)
          {
            continue;
          }

          /**
           * Filter (2): blunt intensity filter
           * Are the mono-isotopic peak intensities of all peptides above the cutoff?
           */
          bool bluntVeto = monoIsotopicPeakIntensityFilter_(patterns_[pattern], spectrum, mz_shifts_actual_indices);
          if (bluntVeto)
          {
            continue;
          }

          /**
           * Filter (3): non-local intensity filter
           * Are the peak intensities of all peptides above the cutoff?
           */
          std::vector<double> intensities_actual; // peak intensities @ m/z peak position + actual m/z shift
          int peaks_found_in_all_peptides_centroided = nonLocalIntensityFilter_(patterns_[pattern], spectrum, mz_shifts_actual_indices, intensities_actual, peaks_found_in_all_peptides);
          if (peaks_found_in_all_peptides_centroided < peaks_per_peptide_min_)
          {
            continue;
          }

          /**
           * Filter (4): zeroth peak filter
           * There should not be a significant peak to the left of the mono-isotopic
           * (i.e. first) peak.
           */
          bool zero_peak = zerothPeakFilter_(patterns_[pattern], intensities_actual);
          if (zero_peak)
          {
            continue;
          }

          /**
           * Filter (5): peptide similarity filter
           * How similar are the isotope patterns of the peptides?
           */
          bool peptide_similarity = peptideSimilarityFilter_(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_centroided);
          if (!peptide_similarity)
          {
            continue;
          }

          /**
           * Filter (6): averagine similarity filter
           * Does each individual isotope pattern resemble a peptide?
           */
          bool averagine_similarity = averagineSimilarityFilter_(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_centroided, peak_position[peak]);
          if (!averagine_similarity)
          {
            continue;
          }

          /**
           * All filters passed.
           */
          // add the peak to the result
          vector<MultiplexFilterResultRaw> results_raw;
          result.addFilterResultPeak(peak_position[peak], rt_picked, mz_shifts_actual, intensities_actual, results_raw);

          // blacklist peaks in the current spectrum and the two neighbouring ones
          blacklistPeaks_(patterns_[pattern], spectrum, mz_shifts_actual_indices, peaks_found_in_all_peptides_centroided);

        }
      }

      // add results of this pattern to list
      filter_results.push_back(result);
    }

    endProgress();

    return filter_results;
  }

  int MultiplexFilteringCentroided::nonLocalIntensityFilter_(const MultiplexIsotopicPeakPattern& pattern, int spectrum_index, const std::vector<int>& mz_shifts_actual_indices, std::vector<double>& intensities_actual, int peaks_found_in_all_peptides) const
  {
    PeakMap::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;

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
}
