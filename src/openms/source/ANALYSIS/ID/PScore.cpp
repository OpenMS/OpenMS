// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/ID/PScore.h>
#include <OpenMS/ANALYSIS/ID/AScore.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>

using std::map;
using std::vector;

namespace OpenMS
{
  vector<Size> PScore::calculateIntensityRankInMZWindow(const vector<double>& mz, const vector<double>& intensities, double mz_window = 100)
  {
    vector<Size> ranks; // note: ranks are zero based
    if (mz.empty())
    {
      return ranks;
    }
    ranks.reserve(mz.size());

    const double half_window = mz_window / 2.0;
    for (Size p = 0; p < mz.size(); ++p)
    {
      const double m = mz[p];
      const double i = intensities[p];

      Size rank(0);

      // count neighbors to the left that have higher intensity
      for (Int j = p - 1; j >= 0; --j)
      {
        if (mz[j] < m - half_window) break;
        if (intensities[j] > i) ++rank;
      }

      // count neighbors to the right that have higher intensity
      for (Size j = p + 1; j < mz.size(); j++)
      {
        if (mz[j] > m + half_window) break;
        if (intensities[j] > i) ++rank;
      }
      ranks.push_back(rank);
    }

    return ranks;
  }


  vector<vector<Size> > PScore::calculateRankMap(const PeakMap& peak_map, double mz_window)
  {
    vector<std::vector<Size> > rank_map; // note: ranks are zero based
    rank_map.reserve(peak_map.size());
    for (Size i = 0; i != peak_map.size(); ++i)
    {
      const PeakSpectrum& spec = peak_map[i];
      vector<double> mz;
      vector<double> intensities;
      for (Size j = 0; j != spec.size(); ++j)
      {
        mz.push_back(spec[j].getMZ());
        intensities.push_back(spec[j].getIntensity());
      }
      rank_map.push_back(calculateIntensityRankInMZWindow(mz, intensities, mz_window));
    }
    return rank_map;
  }

  map<Size, PeakSpectrum > PScore::calculatePeakLevelSpectra(const PeakSpectrum& spec, const vector<Size>& ranks, Size min_level, Size max_level)
  {
    map<Size, MSSpectrum > peak_level_spectra;

    if (spec.empty()) return peak_level_spectra;

    // loop over all peaks and associated (zero-based) ranks
    for (Size i = 0; i != ranks.size(); ++i)
    {
      // start at the highest (less restrictive) level
      for (int j = static_cast<int>(max_level); j >= static_cast<int>(min_level); --j)
      {
        // if the current peak is annotated to have lower or equal rank then allowed for this peak level add it
        if (static_cast<int>(ranks[i]) <= j)
        {
          peak_level_spectra[j].push_back(spec[i]);
        }
        else
        {
          // if the current peak has higher rank than the current level then all it is also to high for the lower levels
          break;
        }
      }
    }
    return peak_level_spectra;
  }

  double PScore::computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const map<Size, PeakSpectrum>& peak_level_spectra, const vector<PeakSpectrum> & theo_spectra, double mz_window)
  {
    AScore a_score_algorithm; // TODO: make the cumulative score function static

    double best_pscore = 0.0;

    for (vector<PeakSpectrum>::const_iterator theo_spectra_it = theo_spectra.begin(); theo_spectra_it != theo_spectra.end(); ++theo_spectra_it)
    {
      const PeakSpectrum& theo_spectrum = *theo_spectra_it;

      // number of theoretical ions for current spectrum
      Size N = theo_spectrum.size();

      for (map<Size, PeakSpectrum>::const_iterator l_it = peak_level_spectra.begin(); l_it != peak_level_spectra.end(); ++l_it)
      {
        const double level = static_cast<double>(l_it->first);
        const PeakSpectrum& exp_spectrum = l_it->second;

        Size matched_peaks(0);
        for (PeakSpectrum::ConstIterator theo_peak_it = theo_spectrum.begin(); theo_peak_it != theo_spectrum.end(); ++theo_peak_it)
        {
          const double& theo_mz = theo_peak_it->getMZ();

          double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);
          double exp_mz = exp_spectrum[index].getMZ();

          // found peak match
          if (std::abs(theo_mz - exp_mz) < max_dist_dalton)
          {
            ++matched_peaks;
          }
        }

        // compute p score as e.g. in the AScore implementation or Andromeda
        const double p = level / mz_window;
        const double pscore = -10.0 * log10(a_score_algorithm.computeCumulativeScore_(N, matched_peaks, p));
        if (pscore > best_pscore)
        {
          best_pscore = pscore;
        }
      }
    }

    return best_pscore;
  }

  double PScore::computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const map<Size, PeakSpectrum>& peak_level_spectra, const PeakSpectrum & theo_spectrum, double mz_window)
  {
    AScore a_score_algorithm; // TODO: make the cumulative score function static

    double best_pscore = 0.0;

    // number of theoretical ions for current spectrum
    Size N = theo_spectrum.size();

    for (map<Size, PeakSpectrum>::const_iterator l_it = peak_level_spectra.begin(); l_it != peak_level_spectra.end(); ++l_it)
    {
      const double level = static_cast<double>(l_it->first);
      const PeakSpectrum& exp_spectrum = l_it->second;

      Size matched_peaks(0);
      if (fragment_mass_tolerance_unit_ppm)
      {
        MatchedIterator<PeakSpectrum, PpmTrait> it(theo_spectrum, exp_spectrum, fragment_mass_tolerance);
        for (; it != it.end(); ++it) ++matched_peaks;
      }
      else
      {
        MatchedIterator<PeakSpectrum, DaTrait> it(theo_spectrum, exp_spectrum, fragment_mass_tolerance);
        for (; it != it.end(); ++it) ++matched_peaks;
      }

      // compute p score as e.g. in the AScore implementation or Andromeda
      const double p = (level + 1) / mz_window;

      const double pscore = -10.0 * log10(a_score_algorithm.computeCumulativeScore_(N, matched_peaks, p));

      if (pscore > best_pscore)
      {
        best_pscore = pscore;
      }
    }

    return best_pscore;
  }

   double massCorrectionTerm(double mass)
   {
     return 0.024 * (mass - 600.0);
   }

   double cleavageCorrectionTerm(Size cleavages, bool consecutive_cleavage)
   {
     switch (cleavages)
     {
       case 0: return 53.2;
       case 1: return consecutive_cleavage ? 42.1 : 31.1;
       case 2: return 17.0;
       default: return 0.0;
     }
   }

   double modificationCorrectionTerm(Size modifications)
   {
     switch (modifications)
     {
       case 0:
         return 42.0;
       case 1:
         return 28.0;
       case 2:
         return 22.0;
       case 3:
         return 16.0;
       case 4:
         return 9.0;
       case 5:
         return 5.0;
       case 6:
         return 2.0;
       default:
         return 0.0;
     }
   }

}

