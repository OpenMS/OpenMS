// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------


#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>

namespace OpenMS
{

// static
void Deisotoper::deisotopeWithAveragineModel(MSSpectrum& spec,
  double fragment_tolerance,
  bool fragment_unit_ppm,
  int number_of_final_peaks,
  int min_charge,
  int max_charge,
  bool keep_only_deisotoped,
  unsigned int min_isopeaks,
  unsigned int max_isopeaks,
  bool make_single_charged,
  bool annotate_charge,
  bool annotate_iso_peak_count,
  bool add_up_intensity)
{
  OPENMS_PRECONDITION(spec.isSorted(), "Spectrum must be sorted.");

    if ((fragment_unit_ppm && fragment_tolerance > 100) || (!fragment_unit_ppm && fragment_tolerance > 0.1))
    {
        throw Exception::IllegalArgument(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "Fragment tolerance must not be greater than 100 ppm or 0.1 Da");
    }

  if (min_isopeaks < 2 || max_isopeaks < 2 || min_isopeaks > max_isopeaks)
  {
    throw Exception::IllegalArgument(__FILE__,
      __LINE__,
      OPENMS_PRETTY_FUNCTION,
      "Minimum/maximum number of isotopic peaks must be at least 2 (and min_isopeaks <= max_isopeaks).");
  }

  if (spec.empty()) { return; }

  // remove 0 intensity peaks
  ThresholdMower threshold_mower_filter;
  threshold_mower_filter.filterPeakSpectrum(spec);

  // discard low-intensity peaks
  if (number_of_final_peaks > 0)
  {
    // only keep number_of_final_peaks highest peaks
    NLargest nlargest_filter = NLargest(number_of_final_peaks);
    nlargest_filter.filterPeakSpectrum(spec);
    
    spec.sortByPosition();
  }

  // store additional information if specified
  Size charge_index(0);
  Size iso_peak_count_index(0);

  // reserve integer data array to store charge of peaks
  if (annotate_charge)
  {
    // expand to hold one additional integer data array to hold the charge
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("charge");
    charge_index = spec.getIntegerDataArrays().size() - 1;
  }
  // reserve integer data array to store number of isotopic peaks for each isotopic pattern
  if (annotate_iso_peak_count)
  {
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("iso_peak_count");
    iso_peak_count_index = spec.getIntegerDataArrays().size() - 1;
  }

  // during discovery phase, work on a constant reference (just to make sure we do not modify spec)
  const MSSpectrum& old_spectrum = spec;

  // determine charge seeds and extend them
  std::vector<Size> mono_isotopic_peak_charge(old_spectrum.size(), 0);
  std::vector<Int> features(old_spectrum.size(), -1);
  std::vector<double> mono_iso_peak_intensity(old_spectrum.size(), 0);
  std::vector<Size> iso_peak_count(old_spectrum.size(), 1);
  int feature_number = 0;

  std::vector<Size> extensions;

  // this should be a vector of vectors of length <= max_isopeaks to reflect multiple possible extensions, but for runtime, this is implemented as a flat vector
  std::vector<Size> clusters;
  // keeps track of the stored clusters sizes
  std::vector<Size> stored_cluster_size;
  std::vector<int> charges_of_extensions;

  const float averagine_check_threshold[7] = {0.0f, 0.0f, 0.05f, 0.1f, 0.2f, 0.4f, 0.6f};

  bool has_precursor_data(false);
  double precursor_mass(0);
  if (old_spectrum.getPrecursors().size() == 1)
  {
    has_precursor_data = true;
    int precursor_charge = old_spectrum.getPrecursors()[0].getCharge();
    precursor_mass = (old_spectrum.getPrecursors()[0].getMZ() * precursor_charge) - (Constants::PROTON_MASS_U * precursor_charge);
  }

  for (size_t current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
  {
    // only process peaks which are not already in a cluster. Would form clusters identical to the cluster they are assigned to,
    // excluding its first peak(s), because only peaks with higher mz are considered for cluster-formation.
    if (features[current_peak] != -1) continue;

    const double current_mz = old_spectrum[current_peak].getMZ();
    const double current_intensity = old_spectrum[current_peak].getIntensity();
    const double tolerance_dalton = fragment_unit_ppm ? Math::ppmToMass(fragment_tolerance, current_mz) : fragment_tolerance;

    clusters.clear();
    stored_cluster_size.clear();
    charges_of_extensions.clear();

    for (int q = min_charge; q <= max_charge; ++q) // all charges are always considered -> order does not matter
    {      
      // try to extend isotopes from mono-isotopic peak

      // do not bother testing charges q (and masses m) with: m/q > precursor_mass/q (or m > precursor_mass)
      if (has_precursor_data)
      {
        double current_theo_mass = (current_mz * q) - (Constants::PROTON_MASS_U * q);
        if (current_theo_mass > (precursor_mass + tolerance_dalton))
        {
          continue;
        }
      }

      // fail early if you do not find a single extension
      const double expected_first_mz = current_mz + Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
      Int p = old_spectrum.findNearest(expected_first_mz, tolerance_dalton);
      if (p == -1) // test for missing peak
      { 
        continue;
      }

      bool found_first_extension = true;

      bool has_min_isopeaks = true;

      extensions.clear();
      extensions.push_back(current_peak);

      // Save frequently used values for performance reasons
      std::vector<double> extensions_intensities = {current_intensity};
      
      // generate averagine distribution for peptide mass corresponding to current mz and charge
      std::vector<double> distr = CoarseIsotopePatternGenerator::approximateIntensities(q * (current_mz - Constants::PROTON_MASS_U), max_isopeaks);
      
      // sum of intensities of both observed and generated peaks is needed for normalization
      double spec_total_intensity = current_intensity;
      double dist_total_intensity = distr[0];

      for (UInt i = 1; i < max_isopeaks; ++i)
      {
        // find next extension (if not first; we found that already)
        if (found_first_extension)
        {
          found_first_extension = false;
        }
        else
        {
          const double expected_mz = current_mz + static_cast<double>(i) * Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
          p = old_spectrum.findNearest(expected_mz, tolerance_dalton);

          if (p == -1)// test for missing peak
          {
            has_min_isopeaks = (i >= min_isopeaks);
            break;
          }
        }

        extensions_intensities.push_back(old_spectrum[p].getIntensity());

        // compare to averagine distribution

        // update sums of intensities
        spec_total_intensity += extensions_intensities.back();
        dist_total_intensity += distr[extensions.size()];

        // compute KL divergence (Sum over all x: P(x) * log(P(x) / Q(x));
        float KL = 0;
        for (unsigned int peak = 0; peak != extensions.size() + 1; ++peak)
        {
          // normalize spectrum intensities and averagine distribution intensities as this is a density measure and 
          // thresholds were probably determined for distributions adding up to 1
          double Px = extensions_intensities[peak] / spec_total_intensity;
          KL += Px * log(Px / (distr[peak] / dist_total_intensity));
        }

        // choose threshold corresponding to cluster size
        float curr_threshold = (extensions.size() + 1 >= 6) ? averagine_check_threshold[6] : averagine_check_threshold[extensions.size() + 1];
          
        // compare to threshold and stop extension if distribution does not fit well enough
        if (KL > curr_threshold)
        {
          has_min_isopeaks = (i >= min_isopeaks);
          break;
        }

        // if model check passed:
        extensions.push_back(p);

        if (annotate_iso_peak_count)
        {
          // with "+ 1" the monoisotopic peak is counted as well
          // current peak is the monoisotopic peak
          iso_peak_count[current_peak] = i + 1;
        }
      } // cluster has been formed

      if (has_min_isopeaks)
      {
        for (UInt i = 0; i != max_isopeaks; ++i)
        {
          clusters.push_back(i < extensions.size() ? extensions[i] : 0);
        }
        stored_cluster_size.push_back(extensions.size());
        charges_of_extensions.push_back(q);
      }
    } // all charges tested, clusters complete
    // if current_peak is possible monoisotopic peak, pick the best of its clusters and annotate peaks with a feature number
    if (!clusters.empty())
    {
      // pick cluster with largest size and highest charge (since all have the same monoisotopic peak)
      UInt best_idx = 0;
      if (stored_cluster_size.size() > 1)// more than one cluster was found
      {
        Size largest_size = 0;
        Int best_cluster_charge = min_charge;

        for (UInt i = 0; i != stored_cluster_size.size(); ++i)
        {
          if ((stored_cluster_size[i] > largest_size) || ((stored_cluster_size[i] == largest_size) && (charges_of_extensions[i] > best_cluster_charge)))
          {
            largest_size = stored_cluster_size[i];
            best_cluster_charge = charges_of_extensions[i];
            best_idx = i;
          }
        }
      }

      // save result
      mono_isotopic_peak_charge[current_peak] = charges_of_extensions[best_idx];
      for (UInt i = 0; i != stored_cluster_size[best_idx]; ++i)
      {
        features[clusters[best_idx * max_isopeaks + i]] = feature_number;

        if (add_up_intensity)
        {
          mono_iso_peak_intensity[current_peak] += old_spectrum[clusters[best_idx * max_isopeaks + i]].getIntensity();
        }
      }
      ++feature_number;
    }
  }

  // apply changes, i.e. select the indices which should survive
  std::vector<Size> select_idx;

  for (size_t i = 0; i != spec.size(); ++i)
  {

    Size z = mono_isotopic_peak_charge[i];
    if (annotate_charge)
    {
      spec.getIntegerDataArrays()[charge_index].push_back((int)z);
    }
    if (annotate_iso_peak_count)
    {
      spec.getIntegerDataArrays()[iso_peak_count_index].push_back((int)iso_peak_count[i]);
    }

    if (!keep_only_deisotoped)
    { // keep all unassigned peaks
      if (features[i] < 0) // this peak is not part of a cluster
      {
        select_idx.push_back(i);
        continue;
      }
    }

    if (z != 0) // this is a monoisotopic peak
    {
      if (add_up_intensity)
      {
        spec[i].setIntensity(mono_iso_peak_intensity[i]);
      }

      // convert mono-isotopic peak with charge assigned by deisotoping
      if (make_single_charged)
      {
        spec[i].setMZ(spec[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
      }
      select_idx.push_back(i);
    }
  }

  // properly subsets all datapoints (incl. dataArrays)
  spec.select(select_idx);
  if (make_single_charged)
  {
    spec.sortByPosition();
  }
  return;
}

// static
void Deisotoper::deisotopeAndSingleCharge(MSSpectrum& spec,
                      double fragment_tolerance,
                      bool fragment_unit_ppm,
                      int min_charge,
                      int max_charge,
                      bool keep_only_deisotoped,
                      unsigned int min_isopeaks,
                      unsigned int max_isopeaks,
                      bool make_single_charged,
                      bool annotate_charge,
                      bool annotate_iso_peak_count,
                      bool use_decreasing_model,
                      unsigned int start_intensity_check,
                      bool add_up_intensity,
                      bool annotate_features)
{
  OPENMS_PRECONDITION(spec.isSorted(), "Spectrum must be sorted.");

    if ((fragment_unit_ppm && fragment_tolerance > 100) || (!fragment_unit_ppm && fragment_tolerance > 0.1))
    {
        throw Exception::IllegalArgument(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "Fragment tolerance must not be greater than 100 ppm or 0.1 Da");
    }

  if (min_isopeaks < 2 || max_isopeaks < 2 || min_isopeaks > max_isopeaks)
  {
    throw Exception::IllegalArgument(__FILE__,
		    __LINE__,
		    OPENMS_PRETTY_FUNCTION,
		    "Minimum/maximum number of isotopic peaks must be at least 2 (and min_isopeaks <= max_isopeaks).");
  }

  if (spec.empty())
  { 
    return; 
  }

  Size charge_index{};
  Size iso_peak_count_index{}, feature_number_dataarray_index{};

  // reserve integer data array to store charge of peaks
  if (annotate_charge)
  {
    // expand to hold one additional integer data array to hold the charge
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("charge");
    charge_index = spec.getIntegerDataArrays().size()-1;
  }
  // reserve integer data array to store number of isotopic peaks for each isotopic pattern
  if (annotate_iso_peak_count)
  {
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("iso_peak_count");
    iso_peak_count_index = spec.getIntegerDataArrays().size()-1;
  }

  if (annotate_features)
  {
    spec.getIntegerDataArrays().resize(spec.getIntegerDataArrays().size() + 1);
    spec.getIntegerDataArrays().back().setName("feature_number");
    feature_number_dataarray_index = spec.getIntegerDataArrays().size() - 1;  
  }

  // during discovery phase, work on a constant reference (just to make sure we do not modify spec)
  const MSSpectrum& old_spectrum = spec;

  // determine charge seeds and extend them
  std::vector<size_t> mono_isotopic_peak(old_spectrum.size(), 0);
  std::vector<int> features(old_spectrum.size(), -1);
  std::vector<double> mono_iso_peak_intensity(old_spectrum.size(), 0);
  std::vector<Size> iso_peak_count(old_spectrum.size(), 1);
  int feature_number = 0;

  std::vector<size_t> extensions;

  bool has_precursor_data(false);
  double precursor_mass(0);
  if (old_spectrum.getPrecursors().size() == 1)
  {
    has_precursor_data = true;
    int precursor_charge = old_spectrum.getPrecursors()[0].getCharge();
    precursor_mass = (old_spectrum.getPrecursors()[0].getMZ() * precursor_charge) - (Constants::PROTON_MASS * precursor_charge);
  }

  for (size_t current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
  {
    const double current_mz = old_spectrum[current_peak].getMZ();
    if (add_up_intensity)
    {
      mono_iso_peak_intensity[current_peak] = old_spectrum[current_peak].getIntensity();
    }

    for (int q = max_charge; q >= min_charge; --q) // important: test charge hypothesis from high to low
    {
      // try to extend isotopes from mono-isotopic peak
      // if extension larger then min_isopeaks possible:
      //   - save charge q in mono_isotopic_peak[]
      //   - annotate_charge all isotopic peaks with feature number
      if (features[current_peak] == -1) // only process peaks which have no assigned feature number
      {
        bool has_min_isopeaks = true;
        const double tolerance_dalton = fragment_unit_ppm ? Math::ppmToMass(fragment_tolerance, current_mz) : fragment_tolerance;

        // do not bother testing charges q (and masses m) with: m/q > precursor_mass/q (or m > precursor_mass)
        if (has_precursor_data)
        {
          double current_theo_mass = (current_mz * q) - (Constants::PROTON_MASS * q);
          if (current_theo_mass > (precursor_mass + tolerance_dalton))
          {
            continue;
          }
        }

        extensions.clear();
        extensions.push_back(current_peak);
        for (unsigned int i = 1; i < max_isopeaks; ++i)
        {
          const double expected_mz = current_mz + static_cast<double>(i) * Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
          const int p = old_spectrum.findNearest(expected_mz, tolerance_dalton);
          if (p == -1) // test for missing peak
          {
            has_min_isopeaks = (i >= min_isopeaks);
            break;
          }
          else
          {
            // Possible improvement: include proper averagine model filtering
            // for now start at the peak with i = start_intensity_check to test hypothesis
            // if start_intensity_check = 0 or 1, start checking by comparing monoisotopic and second isotopic peak
            // if start_intensity_check = 2, start checking by comparing second isotopic peak with the third, etc.
            // Note: this is a common approach used in several other search engines
            if (use_decreasing_model && (i >= start_intensity_check) && (old_spectrum[p].getIntensity() > old_spectrum[extensions.back()].getIntensity()))
            {
              has_min_isopeaks = (i >= min_isopeaks);
              break;
            }
            // averagine check passed or skipped
            extensions.push_back(p);
            if (annotate_iso_peak_count)
            {
              iso_peak_count[current_peak] = i + 1; // with "+ 1" the monoisotopic peak is counted as well
            }
          }
        }

        if (has_min_isopeaks)
        {
          mono_isotopic_peak[current_peak] = q;
          for (unsigned int i = 0; i != extensions.size(); ++i)
          {
            features[extensions[i]] = feature_number;
            // monoisotopic peak intensity is already set above, add up the other intensities here
            if (add_up_intensity && (i != 0))
            {
              mono_iso_peak_intensity[current_peak] += old_spectrum[extensions[i]].getIntensity();
            }
          }
          ++feature_number;
        }
      }
    }
  }

  if (annotate_features)
  { // assign feature indices without copy
    spec.getIntegerDataArrays()[feature_number_dataarray_index].std::vector<Int>::swap(features);
  }
  
  // apply changes, i.e. select the indices which should survive
  std::vector<Size> select_idx;

  for (size_t i = 0; i != spec.size(); ++i)
  {
    Size z = mono_isotopic_peak[i];
    if (annotate_charge)
    {
      spec.getIntegerDataArrays()[charge_index].push_back((int)z);
    }
    if (annotate_iso_peak_count)
    {
      spec.getIntegerDataArrays()[iso_peak_count_index].push_back((int)iso_peak_count[i]);
    }
    if (add_up_intensity)
    {
      spec[i].setIntensity(mono_iso_peak_intensity[i]);
    }

    if (!keep_only_deisotoped)
    { // keep all unassigned peaks
      if (features[i] < 0)
      {
        select_idx.push_back(i);
        continue;
      }
    }

    if (z == 0)
    {
      continue;
    }
    // convert mono-isotopic peak with charge assigned by deisotoping
    if (make_single_charged)
    {
      spec[i].setMZ(spec[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
    }
    select_idx.push_back(i);
  }

  // properly subsets all datapoints (incl. dataArrays)
  spec.select(select_idx);
  spec.sortByPosition();
  return;
}
} // namespace
