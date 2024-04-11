// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include "OpenMS/METADATA/InstrumentSettings.h"
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  
  std::vector<double> PrecursorPurity::computeSingleScanPrecursorPurities(int ms2_spec_idx, int precursor_spec_idx, const MSExperiment & exp, double max_precursor_isotope_deviation)
  {
    const auto& ms2_spec = exp[ms2_spec_idx];
    const auto& precursor_spec = exp[precursor_spec_idx];
    std::vector<double> purities(ms2_spec.getPrecursors().size(), 1.0);
    if (precursor_spec.empty()) return purities; // TODO fail instead?

    Size precursor_idx = 0;
    for (const auto& precursor_info : ms2_spec.getPrecursors())
    {
      typedef PeakMap::SpectrumType::ConstIterator const_spec_iterator;

      // compute distance between isotopic peaks based on the precursor charge.
      const double charge_dist = Constants::NEUTRON_MASS_U / static_cast<double>(precursor_info.getCharge());

      // the actual boundary values
      const double strict_lower_mz = precursor_info.getMZ() - precursor_info.getIsolationWindowLowerOffset();
      const double strict_upper_mz = precursor_info.getMZ() + precursor_info.getIsolationWindowUpperOffset();
      if (strict_lower_mz == strict_upper_mz)
      {
        return purities;
      }

      const double dev_ppm = max_precursor_isotope_deviation / 1e6;
      const double fuzzy_lower_mz = strict_lower_mz * (1 - dev_ppm);
      const double fuzzy_upper_mz = strict_upper_mz * (1 + dev_ppm);

      // first find the actual precursor peak
      Size precursor_peak_idx = precursor_spec.findNearest(precursor_info.getMZ());
      const Peak1D& precursor_peak = precursor_spec[precursor_peak_idx];

      // now we get ourselves some border iterators
      const_spec_iterator lower_bound = precursor_spec.MZBegin(fuzzy_lower_mz);
      const_spec_iterator upper_bound = precursor_spec.MZEnd(precursor_info.getMZ());

      Peak1D::IntensityType precursor_intensity = precursor_peak.getIntensity();
      Peak1D::IntensityType total_intensity = precursor_peak.getIntensity();

      // ------------------------------------------------------------------------------
      // try to find a match for our isotopic peak on the left side

      double expected_next_mz = precursor_peak.getMZ() - charge_dist;

      while (expected_next_mz > fuzzy_lower_mz)
      {
        // find nearest peak in precursor window
        const_spec_iterator np_it = precursor_spec.MZBegin(lower_bound, expected_next_mz, upper_bound);

        // handle border cases

        // check if next peak has smaller dist
        const_spec_iterator np_it2 = np_it;
        ++np_it;

        if (std::fabs(np_it2->getMZ() - expected_next_mz) < std::fabs(np_it->getMZ() - expected_next_mz))
        {
          np_it = np_it2;
        }

        // compute difference between found peak and expected
        double min_ppm_diff = std::fabs(np_it->getMZ() - expected_next_mz)  * 1000000 / expected_next_mz;

        // check if we found an isotopic peak
        if (min_ppm_diff < max_precursor_isotope_deviation)
        {
          if (np_it->getMZ() > strict_lower_mz)
          {
            precursor_intensity += np_it->getIntensity();
          }
          else
          {
            // we're in the fuzzy area, so we will take only 50% of the given intensity
            // since we assume that the isolation window borders are not sharp
            precursor_intensity += 0.5 * np_it->getIntensity();
          }

          // update expected_next_mz
          expected_next_mz = np_it->getMZ() - charge_dist;
        }
        else
        {
          // update expected_next_mz with theoretical position
          expected_next_mz -= charge_dist;
        }
      }

      // ------------------------------------------------------------------------------
      // try to find a match for our isotopic peak on the right

      // redefine bounds
      lower_bound = precursor_spec.MZBegin(precursor_info.getMZ());
      upper_bound = precursor_spec.MZEnd(fuzzy_upper_mz);

      expected_next_mz = precursor_peak.getMZ() + charge_dist;

      while (expected_next_mz < fuzzy_upper_mz)
      {
        // find nearest peak in precursor window
        const_spec_iterator np_it = precursor_spec.MZBegin(lower_bound, expected_next_mz, upper_bound);

        // handle border cases

        // check if next peak has smaller dist
        const_spec_iterator np_it2 = np_it;
        ++np_it;

        if (std::fabs(np_it2->getMZ() - expected_next_mz) < std::fabs(np_it->getMZ() - expected_next_mz))
        {
          np_it = np_it2;
        }

        // compute difference between found peak and expected
        double min_ppm_diff = std::fabs(np_it->getMZ() - expected_next_mz)  * 1000000 / expected_next_mz;

        // check if we found an isotopic peak
        if (min_ppm_diff < max_precursor_isotope_deviation)
        {
          if (np_it->getMZ() < strict_upper_mz)
          {
            precursor_intensity += np_it->getIntensity();
          }
          else
          {
            // we're in the fuzzy area, so we will take only 50% of the given intensity
            // since we assume that the isolation window borders are not sharp
            precursor_intensity += 0.5 * np_it->getIntensity();
          }

          // update expected_next_mz
          expected_next_mz = np_it->getMZ() + charge_dist;
        }
        else
        {
          // update expected_next_mz with theoretical position
          expected_next_mz += charge_dist;
        }
      }

      // ------------------------------------------------------------------------------
      // compute total intensity
      int idx = static_cast<int>(precursor_peak_idx) - 1;
      while (idx >= 0 && precursor_spec[idx].getMZ() > fuzzy_lower_mz)
      {
        if (precursor_spec[idx].getMZ() > strict_lower_mz)
        {
          total_intensity += precursor_spec[idx].getIntensity();
        }
        else
        {
          // we're in the fuzzy area, so we will take only 50% of the given intensity
          // since we assume that the isolation window borders are not sharp
          total_intensity += 0.5 * precursor_spec[idx].getIntensity();
        }
        --idx;
      }

      idx = static_cast<int>(precursor_peak_idx) + 1;
      while (idx < static_cast<int>(precursor_spec.size()) && precursor_spec[idx].getMZ() < fuzzy_upper_mz)
      {
        if (precursor_spec[idx].getMZ() < strict_upper_mz)
        {
          total_intensity += precursor_spec[idx].getIntensity();
        }
        else
        {
          // we're in the fuzzy area, so we will take only 50% of the given intensity
          // since we assume that the isolation window borders are not sharp
          total_intensity += 0.5 * precursor_spec[idx].getIntensity();
        }
        ++idx;
      }

      purities[precursor_idx] = precursor_intensity / total_intensity;
      precursor_idx++;
    }
    return purities;
  }

  std::vector<double> PrecursorPurity::computeInterpolatedPrecursorPurity(int ms2_spec_idx, int precursor_spec_idx, int next_ms1_spec_idx, const MSExperiment & exp, double max_precursor_isotope_deviation)
  {
    const auto& ms2_spec = exp[ms2_spec_idx];
    const auto& precursor_spec = exp[precursor_spec_idx];
    const auto& next_ms1_spec = exp[next_ms1_spec_idx];
    // compute purity of preceding ms1 scan
    std::vector<double> early_scan_purity = computeSingleScanPrecursorPurities(ms2_spec_idx, precursor_spec_idx, exp, max_precursor_isotope_deviation);
    std::vector<double> late_scan_purity  = computeSingleScanPrecursorPurities(ms2_spec_idx, next_ms1_spec_idx, exp, max_precursor_isotope_deviation);
    std::vector<double> interpolated_purity;
    interpolated_purity.reserve(early_scan_purity.size());
    for (Size i = 0; i < early_scan_purity.size(); ++i)
    {
      // calculating the extrapolated, S2I value as a time weighted linear combination of the two scans
      // see: Savitski MM, Sweetman G, Askenazi M, Marto JA, Lang M, Zinn N, et al. (2011).
      // Analytical chemistry 83: 8959–67. http://www.ncbi.nlm.nih.gov/pubmed/22017476
      // std::fabs is applied to compensate for potentially negative RTs
      interpolated_purity.push_back(
        std::fabs(ms2_spec.getRT() - precursor_spec.getRT()) *
            ((late_scan_purity[i] - early_scan_purity[i]) / std::fabs(next_ms1_spec.getRT() - precursor_spec.getRT()))
            + early_scan_purity[i]);
    }
    return interpolated_purity;
  }

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
