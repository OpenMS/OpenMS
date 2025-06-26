// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/NUXL/NuXLPreprocessSpectra.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLDeisotoper.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLConstants.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <unordered_set>

namespace OpenMS
{

// static
void NuXLPreprocessSpectra::preprocessSpectra(MSExperiment& exp, 
                                               bool single_charge_spectra, 
                                               bool annotate_charge,
                                               double window_size,
                                               Size peakcount,
                                               const std::map<String, PrecursorPurity::PurityScores>& purities)
{
  // filter MS2 map
  // remove 0 intensities
  ThresholdMower threshold_mower_filter;
  threshold_mower_filter.filterPeakMap(exp);

#pragma omp parallel for
  for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
  {
    MSSpectrum & spec = exp[exp_index];

    // sort by mz
    spec.sortByPosition();

    // deisotope
    NuXLDeisotoper::deisotopeAndSingleCharge(spec, 
                                       0.01,
                                       false,
                                       1, 3, 
                                       false, 
                                       2, 10, 
                                       single_charge_spectra, 
                                       annotate_charge,
                                       false, // no iso peak count annotation
                                       true, // decreasing isotope model
                                       2, // enforce only starting from second peak
                                       true); // add up intensities
  }

  filterPeakInterference_(exp, purities);

  // remove empty spectra as they can cause trouble downstream
  auto& sp = exp.getSpectra();
  sp.erase(std::remove_if(sp.begin(), sp.end(), [](const MSSpectrum& s) { return s.empty(); }), exp.end());

  Normalizer normalizer;
  normalizer.filterPeakMap(exp);

  // sort by rt
  exp.sortSpectra(false);

  // filter settings
  WindowMower window_mower_filter;
  Param filter_param = window_mower_filter.getParameters();
  filter_param.setValue("windowsize", window_size, "The size of the sliding window along the m/z axis.");
  filter_param.setValue("peakcount", peakcount, "The number of peaks that should be kept.");
  filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
  window_mower_filter.setParameters(filter_param);

  NLargest nlargest_filter = NLargest(400);

#pragma omp parallel for
  for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
  {
    MSSpectrum & spec = exp[exp_index];
    // sort by mz
    spec.sortByPosition();

    if (annotate_charge)
    { 
      // set Unknown charge to z=1. Otherwise we get a lot of spurious matches 
      // to highly charged fragments in the low m/z region
      if (spec.empty()) continue; // TODO: maybe add empty integerdataarray in deisotoping? seems to be missing here
      DataArrays::IntegerDataArray& ia = spec.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX]; // charge array
      for (int & z : ia) { if (z == 0) { z = 1; } }
    } 
  
    // remove noise
    window_mower_filter.filterPeakSpectrum(spec);
    
    nlargest_filter.filterPeakSpectrum(spec);
 
    // sort (nlargest changes order)
    spec.sortByPosition();

    // calculate TIC and store in float data array
    spec.getFloatDataArrays().clear();
    spec.getFloatDataArrays().resize(1);
    double TIC = spec.calculateTIC();
    spec.getFloatDataArrays()[0].push_back(TIC);
    spec.getFloatDataArrays()[0].setName("TIC");
  }
}

// static
void NuXLPreprocessSpectra::filterPeakInterference_(MSExperiment& spectra, 
                                                     const std::map<String, PrecursorPurity::PurityScores>& purities, 
                                                     double fragment_mass_tolerance, 
                                                     bool fragment_mass_tolerance_unit_ppm)
{
  double filtered_peaks_count{0};
  size_t filtered_spectra{0};
  for (auto& s : spectra)
  {
    std::unordered_set<size_t> idx_to_remove;
    auto it = purities.find(s.getNativeID());
    if (it != purities.end())
    {        
      for (const auto& interfering_peak : it->second.interfering_peaks)
      {
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? interfering_peak.getMZ()  * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
        auto pos = s.findNearest(interfering_peak.getMZ(), max_dist_dalton, max_dist_dalton); 
        if (pos != -1) 
        {
          idx_to_remove.insert(pos);
        }
      }
      std::vector<size_t> idx_to_keep; // inverse
      for (size_t i = 0; i != s.size(); ++i)
      { // add indices we don't want to remove
        if (idx_to_remove.find(i) == idx_to_remove.end()) idx_to_keep.push_back(i);
      }
      filtered_peaks_count += idx_to_remove.size();
      s.select(idx_to_keep);
    }
    ++filtered_spectra;
  } 
  OPENMS_LOG_INFO << "Filtered out " << filtered_peaks_count << " peaks in total that matched to precursor interference." << std::endl;
  if (filtered_spectra > 0) OPENMS_LOG_INFO << "  On average " << filtered_peaks_count / (double)filtered_spectra << " peaks per MS2." << std::endl;
}

}