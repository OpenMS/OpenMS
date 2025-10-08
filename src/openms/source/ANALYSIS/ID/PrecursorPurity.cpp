// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  PrecursorPurity::PurityScores PrecursorPurity::computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm)
  {
    PrecursorPurity::PurityScores score;
    double target_mz = pre.getMZ();
    double lower = target_mz - pre.getIsolationWindowLowerOffset();
    double upper = target_mz + pre.getIsolationWindowUpperOffset();

    int charge = abs(pre.getCharge());

    if (charge == 0) charge = 1; // prevent division by zero

    double precursor_tolerance_abs = precursor_mass_tolerance_unit_ppm ? (target_mz * precursor_mass_tolerance*2 * 1e-6) : precursor_mass_tolerance*2;

    auto lower_it = ms1.MZBegin(lower);
    auto upper_it = ms1.MZEnd(upper);

    PeakSpectrum isolated_window;
    while (lower_it != upper_it)
    {
      isolated_window.push_back(*lower_it);
      lower_it++;
    }

    // total intensity in isolation window
    double total_intensity(0);
    double target_intensity(0);
    Size target_peak_count(0);
    for (const auto& peak : isolated_window)
    {
      total_intensity += peak.getIntensity();
    }

    // search for the target peak, return scores with 0-values if it is not found
    if (isolated_window.empty())
    {
      return score;
    }

    // estimate a lower boundary for isotopic peaks
    int negative_isotopes((pre.getIsolationWindowLowerOffset() * charge));
    double iso = -negative_isotopes;
    // depending on the isolation window, the first estimated peak might be outside the window
    if (target_mz + (iso * Constants::C13C12_MASSDIFF_U / charge) < lower)
    {
      iso++;
    }

    // deisotoping (try to find isotopic peaks of the precursor mass, even if the actual precursor peak is missing)
    while (true) // runs as long as the next mz is within the isolation window
    {
      double next_peak = target_mz + (iso * Constants::C13C12_MASSDIFF_U / charge);

      // stop loop when new mz is outside the isolation window
      // changes through the isotope index iso
      if (next_peak > upper)
      {
        break;
      }
      int next_iso_index = isolated_window.findNearest(next_peak, precursor_tolerance_abs);
      if (next_iso_index != -1)
      {
        target_intensity += isolated_window[next_iso_index].getIntensity();

        isolated_window.erase(isolated_window.begin() + next_iso_index);
        target_peak_count++;
      }
      // always increment iso to progress the loop
      iso++;
    }

    double rel_sig(0);
    if (target_intensity > 0.0)
    {
      rel_sig = target_intensity / total_intensity;
    }

    score.total_intensity = total_intensity;
    score.target_intensity = target_intensity;
    score.signal_proportion = rel_sig;
    score.target_peak_count = target_peak_count;
    score.interfering_peak_count = isolated_window.size();
    score.interfering_peaks = isolated_window;

    return score;
  }

  PrecursorPurity::PurityScores PrecursorPurity::combinePrecursorPurities(const PrecursorPurity::PurityScores& score1, const PrecursorPurity::PurityScores& score2)
  {
    PrecursorPurity::PurityScores score;
    score.total_intensity = score1.total_intensity + score2.total_intensity;
    score.target_intensity = score1.target_intensity + score2.target_intensity;
    if (score.target_intensity > 0.0) // otherwise default value of 0 is used
    {
      score.signal_proportion = score.target_intensity / score.total_intensity;
    }
    score.target_peak_count = score1.target_peak_count + score2.target_peak_count;
    score.interfering_peak_count = score1.interfering_peak_count + score2.interfering_peak_count;

    return score;
  }

  std::map<String, PrecursorPurity::PurityScores> PrecursorPurity::computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm, bool ignore_missing_precursor_spectra)
  {
    std::map<String, PrecursorPurity::PurityScores> purityscores;
    std::pair<std::map<String, PrecursorPurity::PurityScores>::iterator, bool> insert_return_value;
    int spectra_size = static_cast<int>(spectra.size());

    if (spectra[0].getMSLevel() != 1 && !ignore_missing_precursor_spectra)
    {
      OPENMS_LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. First Spectrum is not MS1. Precursor Purity info will not be calculated!\n";
      return purityscores;
    }

    for (int i = 0; i < spectra_size; ++i)
    {
      if (spectra[i].getMSLevel() == 2)
      {
        auto parent_spectrum_it = spectra.getPrecursorSpectrum(spectra.begin()+i);
        if (parent_spectrum_it == spectra.end() && !ignore_missing_precursor_spectra)
        {
          OPENMS_LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. An MS2 spectrum without parent spectrum detected. Precursor Purity info will not be calculated!\n";
          return std::map<String, PrecursorPurity::PurityScores>();
        }
        if (spectra[i].getNativeID().empty())
        {
          OPENMS_LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. Spectrum without an ID. Precursor Purity info will not be calculated!\n";
          return std::map<String, PrecursorPurity::PurityScores>();
        }

        // check for uniqueness of IDs by inserting initialized (0-value) scores into map
        insert_return_value = purityscores.insert(std::pair<String, PrecursorPurity::PurityScores>(spectra[i].getNativeID(), PrecursorPurity::PurityScores()));
        if (!insert_return_value.second)
        {
          OPENMS_LOG_WARN << "Warning: Input data not suitable for Precursor Purity computation. Duplicate Spectrum IDs. Precursor Purity info will not be calculated!\n";
          return std::map<String, PrecursorPurity::PurityScores>();
        }
      }
    }

#pragma omp parallel for schedule(guided)
    for (int i = 0; i < spectra_size; ++i)
    {
      if (spectra[i].getMSLevel() == 2)
      {
        PrecursorPurity::PurityScores score;
        auto parent_spectrum_it = spectra.getPrecursorSpectrum(spectra.begin() + i);
        if (parent_spectrum_it != spectra.end())
        {
          score = PrecursorPurity::computePrecursorPurity((*parent_spectrum_it), spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);
        }
#pragma omp critical (purityscores_access)
        {
          // replace the initialized values
          purityscores[spectra[i].getNativeID()] = score;
        }
      } // end of MS2 spectrum
    } // end of parallelized spectra loop
    return purityscores;
  } // end of function def

}
