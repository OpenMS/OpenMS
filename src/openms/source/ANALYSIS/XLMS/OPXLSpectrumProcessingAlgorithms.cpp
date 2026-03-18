// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>

// preprocessing and filtering
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace OpenMS
{

  PeakSpectrum OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum)
  {
    // merge peaks: create new spectrum, insert peaks from first and then from second spectrum
    PeakSpectrum resulting_spectrum;
    resulting_spectrum.insert(resulting_spectrum.end(), first_spectrum.begin(), first_spectrum.end());
    resulting_spectrum.insert(resulting_spectrum.end(), second_spectrum.begin(), second_spectrum.end());

    // merge DataArrays in a similar way
    for (Size i = 0; i < first_spectrum.getFloatDataArrays().size(); i++)
    {
      // TODO instead of this "if", get second array by name if available.  would not be dependent on order.
      if (second_spectrum.getFloatDataArrays().size() > i)
      {
        PeakSpectrum::FloatDataArray float_array;
        float_array.insert(float_array.end(), first_spectrum.getFloatDataArrays()[i].begin(), first_spectrum.getFloatDataArrays()[i].end());
        float_array.insert(float_array.end(), second_spectrum.getFloatDataArrays()[i].begin(), second_spectrum.getFloatDataArrays()[i].end());
        resulting_spectrum.getFloatDataArrays().push_back(float_array);
        resulting_spectrum.getFloatDataArrays()[i].setName(first_spectrum.getFloatDataArrays()[i].getName());
      }
    }

    for (Size i = 0; i < first_spectrum.getStringDataArrays().size(); i++)
    {
      if (second_spectrum.getStringDataArrays().size() > i)
      {
        PeakSpectrum::StringDataArray string_array;
        string_array.insert(string_array.end(), first_spectrum.getStringDataArrays()[i].begin(), first_spectrum.getStringDataArrays()[i].end());
        string_array.insert(string_array.end(), second_spectrum.getStringDataArrays()[i].begin(), second_spectrum.getStringDataArrays()[i].end());
        resulting_spectrum.getStringDataArrays().push_back(string_array);
        resulting_spectrum.getStringDataArrays()[i].setName(first_spectrum.getStringDataArrays()[i].getName());
      }
    }

    for (Size i = 0; i < first_spectrum.getIntegerDataArrays().size(); i++)
    {
      if (second_spectrum.getIntegerDataArrays().size() > i)
      {
        PeakSpectrum::IntegerDataArray integer_array;
        integer_array.insert(integer_array.end(), first_spectrum.getIntegerDataArrays()[i].begin(), first_spectrum.getIntegerDataArrays()[i].end());
        integer_array.insert(integer_array.end(), second_spectrum.getIntegerDataArrays()[i].begin(), second_spectrum.getIntegerDataArrays()[i].end());
        resulting_spectrum.getIntegerDataArrays().push_back(integer_array);
        resulting_spectrum.getIntegerDataArrays()[i].setName(first_spectrum.getIntegerDataArrays()[i].getName());
      }
    }

    // Spectra were simply concatenated, so they are not sorted by position anymore
    resulting_spectrum.sortByPosition();
    return resulting_spectrum;
  }

  PeakMap OPXLSpectrumProcessingAlgorithms::preprocessSpectra(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, Size peptide_min_size, Int min_precursor_charge, Int max_precursor_charge, bool deisotope, bool labeled)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);
    OPENMS_LOG_DEBUG << "Deisotoping and filtering spectra." << endl;

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    PeakMap filtered_spectra;


#pragma omp parallel for
    for (SignedSize exp_index = 0; exp_index < static_cast<SignedSize>(exp.size()); ++exp_index)
    {
      // for labeled experiments, the pairs of heavy and light spectra are linked by spectra indices from the consensusXML, so the returned number of spectra has to be equal to the input
      bool process_this_spectrum(labeled);
      if (exp[exp_index].getMSLevel() != 2)
      {
        continue;
      }

      const vector<Precursor>& precursor = exp[exp_index].getPrecursors();

      if (!process_this_spectrum && precursor.size() == 1 && exp[exp_index].size() >= peptide_min_size * 2)
      {
        int precursor_charge = precursor[0].getCharge();
        if (precursor_charge >= min_precursor_charge && precursor_charge <= max_precursor_charge)
        {
          process_this_spectrum = true;
        }
      }

      if (!process_this_spectrum)
      {
        continue;
      }

      if (deisotope)
      {
        PeakSpectrum deisotoped = exp[exp_index];
        Deisotoper::deisotopeAndSingleCharge(deisotoped,
          fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
          1, 7,   // min / max charge
          false,  // keep only deisotoped
          3, 10,  // min / max isopeaks
          false,  // make single charged
          true,   // annotate charge
          true,   // annotate isotopic peak counts
          true,   // use simple averagine model
          3,      // peak to start averagine model
          true    // add up intensity into monoisotopic peak
          );

        // only consider spectra, that have at least as many peaks as two times the minimal peptide size after deisotoping
        if (deisotoped.size() > peptide_min_size * 2 || labeled)
        {
          window_mower_filter.filterPeakSpectrum(deisotoped);
          deisotoped.sortByPosition();

#pragma omp critical (filtered_spectra_access)
          filtered_spectra.addSpectrum(deisotoped);
        }
      }
      else
      {
        PeakSpectrum filtered = exp[exp_index];
        if (!labeled) // this kind of filtering is not necessary for labeled cross-links, since they area filtered by comparing heavy and light spectra later
        {
          window_mower_filter.filterPeakSpectrum(filtered);
        }

        // only consider spectra, that have at least as many peaks as two times the minimal peptide size after filtering
        if (filtered.size() > peptide_min_size * 2 || labeled)
        {
          filtered.sortByPosition();

#pragma omp critical (filtered_spectra_access)
          filtered_spectra.addSpectrum(filtered);
        }
      }
    } // end of parallelized loop over spectra
    return filtered_spectra;
  }

  void OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(
    std::vector<std::pair<Size, Size> > & alignment, double fragment_mass_tolerance,
    bool fragment_mass_tolerance_unit_ppm,
    const PeakSpectrum& theo_spectrum,
    const PeakSpectrum& exp_spectrum,
    const DataArrays::IntegerDataArray& theo_charges,
    const DataArrays::IntegerDataArray& exp_charges,
    DataArrays::FloatDataArray& ppm_error_array,
    double intensity_cutoff)
  {
    OPENMS_PRECONDITION(exp_spectrum.isSorted(), "Spectrum needs to be sorted.");
    OPENMS_PRECONDITION(theo_spectrum.isSorted(), "Spectrum needs to be sorted.");
    OPENMS_PRECONDITION((alignment.empty() == true), "Alignment result vector needs to be empty.");
    OPENMS_PRECONDITION((ppm_error_array.empty() == true), "ppm error result vector needs to be empty.");

    const Size n_t(theo_spectrum.size());
    const Size n_e(exp_spectrum.size());
    const bool has_charge = !(exp_charges.empty() || theo_charges.empty());

    if (n_t == 0 || n_e == 0) { return; }

    Size t(0), e(0);
    alignment.reserve(theo_spectrum.size());
    ppm_error_array.reserve(theo_spectrum.size());

    while (t < n_t && e < n_e)
    {
      const double theo_mz = theo_spectrum[t].getMZ();
      const double exp_mz = exp_spectrum[e].getMZ();

      int tz(0), ez(0);
      if (has_charge)
      {
        tz = theo_charges[t];
        ez = exp_charges[e];
      }
      const bool tz_matches_ez = (ez == tz || !ez || !tz);

      double ti = theo_spectrum[t].getIntensity();
      double ei = exp_spectrum[e].getIntensity();
      const bool initial_intensity_matches = ( std::min(ti, ei) / std::max(ti, ei) ) > intensity_cutoff;


      double d = exp_mz - theo_mz;
      const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

      if (fabs(d) <= max_dist_dalton) // match in tolerance window?
      {
        // get first peak with matching charge in tolerance window
        if (!tz_matches_ez || !initial_intensity_matches)
        {
          Size e_candidate(e);
          while (e_candidate < n_e-1)
          {
            ++e_candidate;
            double new_ez = has_charge ? exp_charges[e_candidate] : 0;
            double new_ei = exp_spectrum[e_candidate].getIntensity();
            const bool charge_matches = (new_ez == tz || !new_ez || !tz);
            const bool intensity_matches = ( std::min(ti, new_ei) / std::max(ti, new_ei) ) > intensity_cutoff;
            double new_d = exp_spectrum[e_candidate].getMZ() - theo_mz;
            if (charge_matches && new_d <= max_dist_dalton && intensity_matches)
            { // found a match
              break;
            }
            else if (new_d > max_dist_dalton)
            { // no match found
              e_candidate = e;
              break;
            }
          }
          if (e == e_candidate)
          { // no match found continue with next theo. peak
            ++t;
            continue;
          }
          else
          { // match found
            e = e_candidate;
          }
        }

        // Invariant: e now points to the first peak in tolerance window, that matches in charge and intensity

        // last peak? there can't be a better one in this tolerance window
        if (e >= n_e - 1)
        {
          // add match
          alignment.emplace_back(std::make_pair(t, e));
          // add ppm error
          double ppm_error = (exp_spectrum[e].getMZ() - theo_mz) / theo_mz * 1e6;
          ppm_error_array.emplace_back(ppm_error);
          return;
        }

        Size closest_exp_peak(e);

        // Invariant: closest_exp_peak always point to best match

        double new_ez(0);
        double best_d = exp_spectrum[closest_exp_peak].getMZ() - theo_mz;

        do // check for better match in tolerance window
        {
          // advance to next exp. peak
          ++e;

          // determine distance of next peak
          double new_d = exp_spectrum[e].getMZ() - theo_mz;
          const bool in_tolerance_window = (fabs(new_d) < max_dist_dalton);

          if (!in_tolerance_window) { break; }

          // Invariant: e is in tolerance window

          // check if charge and intensity of next peak matches
          if (has_charge) { new_ez = exp_charges[e]; }
          const bool charge_matches = (new_ez == tz || !new_ez || !tz);
          double new_ei = exp_spectrum[e].getIntensity();
          const bool intensity_matches = ( std::min(ti, new_ei) / std::max(ti, new_ei) ) > intensity_cutoff;
          if (!charge_matches || !intensity_matches) { continue; }

          // Invariant: charge and intensity matches

          const bool better_distance = (fabs(new_d) <= fabs(best_d));

          // better distance (and matching charge)? better match found
          if (better_distance)
          { // found a better match
            closest_exp_peak = e;
            best_d = new_d;
          }
          else
          { // distance got worse -> no additional matches!
            break;
          }
        }
        while (e < n_e - 1);

        // search in tolerance window for an experimental peak closer to theoretical one
        alignment.emplace_back(std::make_pair(t, closest_exp_peak));

        // add ppm error for this match
        double ppm_error = (exp_spectrum[closest_exp_peak].getMZ() - theo_mz) / theo_mz * 1e6;
        ppm_error_array.emplace_back(ppm_error);

        e = closest_exp_peak + 1;  // advance experimental peak to 1-after the best match
        ++t; // advance theoretical peak
      }
      else if (d < 0) // exp. peak is left of theo. peak (outside of tolerance window)
      {
        ++e;
      }
      else if (d > 0) // theo. peak is left of exp. peak (outside of tolerance window)
      {
        ++t;
      }
    }
  }

  void OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(
    std::vector<std::pair<Size, Size> > & alignment,
    double fragment_mass_tolerance,
    bool fragment_mass_tolerance_unit_ppm,
    const std::vector< SimpleTSGXLMS::SimplePeak >& theo_spectrum,
    const PeakSpectrum& exp_spectrum,
    const DataArrays::IntegerDataArray& exp_charges)
  {
    alignment.clear();
    const Size n_t(theo_spectrum.size());
    const Size n_e(exp_spectrum.size());
    const bool has_charge = !(exp_charges.empty());

    if (n_t == 0 || n_e == 0) { return; }

    Size t(0), e(0);
    alignment.reserve(theo_spectrum.size());

    while (t < n_t && e < n_e)
    {
      const double theo_mz = theo_spectrum[t].mz;
      const double exp_mz = exp_spectrum[e].getMZ();

      int tz(0), ez(0);
      if (has_charge)
      {
        tz = theo_spectrum[t].charge;
        ez = exp_charges[e];
      }
      const bool tz_matches_ez = (ez == tz || !ez || !tz);

      double d = exp_mz - theo_mz;
      const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

      if (fabs(d) <= max_dist_dalton) // match in tolerance window?
      {
        // get first peak with matching charge in tolerance window
        if (!tz_matches_ez)
        {
          Size e_candidate(e);
          while (e_candidate < n_e-1)
          {
            ++e_candidate;
            double new_ez = has_charge ? exp_charges[e_candidate] : 0;
            const bool charge_matches = (new_ez == tz || !new_ez || !tz);
            double new_d = exp_spectrum[e_candidate].getMZ() - theo_mz;
            if (charge_matches && new_d <= max_dist_dalton)
            { // found a match
              break;
            }
            else if (new_d > max_dist_dalton)
            { // no match found
              e_candidate = e;
              break;
            }
          }
          if (e == e_candidate)
          { // no match found continue with next theo. peak
            ++t;
            continue;
          }
          else
          { // match found
            e = e_candidate;
          }
        }

        // Invariant: e now points to the first peak in tolerance window, that matches in charge and intensity

        // last peak? there can't be a better one in this tolerance window
        if (e >= n_e - 1)
        {
          // add match
          alignment.emplace_back(std::make_pair(t, e));
          return;
        }

        Size closest_exp_peak(e);

        // Invariant: closest_exp_peak always point to best match

        double new_ez(0);
        double best_d = exp_spectrum[closest_exp_peak].getMZ() - theo_mz;

        do // check for better match in tolerance window
        {
          // advance to next exp. peak
          ++e;

          // determine distance of next peak
          double new_d = exp_spectrum[e].getMZ() - theo_mz;
          const bool in_tolerance_window = (fabs(new_d) < max_dist_dalton);

          if (!in_tolerance_window) { break; }

          // Invariant: e is in tolerance window

          // check if charge and intensity of next peak matches
          if (has_charge) { new_ez = exp_charges[e]; }
          const bool charge_matches = (new_ez == tz || !new_ez || !tz);
          if (!charge_matches) { continue; }

          // Invariant: charge and intensity matches

          const bool better_distance = (fabs(new_d) <= fabs(best_d));

          // better distance (and matching charge)? better match found
          if (better_distance)
          { // found a better match
            closest_exp_peak = e;
            best_d = new_d;
          }
          else
          { // distance got worse -> no additional matches!
            break;
          }
        }
        while (e < n_e - 1);

        // search in tolerance window for an experimental peak closer to theoretical one
        alignment.emplace_back(std::make_pair(t, closest_exp_peak));
        e = closest_exp_peak + 1;  // advance experimental peak to 1-after the best match
        ++t; // advance theoretical peak
      }
      else if (d < 0) // exp. peak is left of theo. peak (outside of tolerance window)
      {
        ++e;
      }
      else if (d > 0) // theo. peak is left of exp. peak (outside of tolerance window)
      {
        ++t;
      }
    }
  }

}
