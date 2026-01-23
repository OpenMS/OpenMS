// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/SRMHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

std::vector<OpenMS::MSChromatogram> SRMHandler::collectChromatogramsForIrt(const std::vector< OpenSwath::SwathMap > & swath_maps)
{
  std::vector<OpenMS::MSChromatogram> irt_chromatograms;
  size_t total_chroms = 0;
  for (size_t sm_idx = 0; sm_idx < swath_maps.size(); ++sm_idx)
  {
    const auto& sm = swath_maps[sm_idx];
    OpenSwath::SpectrumAccessPtr chrom_map = sm.sptr;
    size_t n_chroms = chrom_map ? chrom_map->getNrChromatograms() : 0;
    for (Size i = 0; i < n_chroms; ++i)
    {
      auto* openms_map = dynamic_cast<OpenMS::SpectrumAccessOpenMS*>(chrom_map.get());
      if (!openms_map) continue;
      auto chrom_ptr = openms_map->getChromatogramById(i);
      if (!chrom_ptr) continue;
      MSChromatogram chrom;
      OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chrom_ptr, chrom);
      String native_id = openms_map->getChromatogramNativeID(i);
      // For iRT calibration, preserve the original nativeID in the chromatogram
      chrom.setNativeID(native_id);
      irt_chromatograms.push_back(chrom);
      ++total_chroms;
    }
  }
  OPENMS_LOG_DEBUG << "SRM/MRM: collected " << total_chroms << " chromatograms for iRT." << std::endl;
  return irt_chromatograms;
}

std::vector<OpenMS::MSChromatogram> SRMHandler::extractAndMapChromatogramsForTransitions(
  const std::vector< OpenSwath::SwathMap > & swath_maps,
  const OpenSwath::LightTargetedExperiment & transition_exp,
  const ChromExtractParams & cp)
{
  // Collect all chromatograms first and store original native ids as meta values
  std::vector<MSChromatogram> all_chroms;
  for (const auto & sm : swath_maps)
  {
    OpenSwath::SpectrumAccessPtr chrom_map = sm.sptr;
    if (!chrom_map) continue;
    auto* openms_map = dynamic_cast<OpenMS::SpectrumAccessOpenMS*>(chrom_map.get());
    if (!openms_map) continue;
    Size n_chroms = openms_map->getNrChromatograms();
    for (Size i = 0; i < n_chroms; ++i)
    {
      auto chrom_ptr = openms_map->getChromatogramById(i);
      if (!chrom_ptr) continue;
      MSChromatogram chrom;
      OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chrom_ptr, chrom);
      String orig = openms_map->getChromatogramNativeID(i);
      chrom.setMetaValue("original_native_id", orig);
      // clear nativeID until mapped
      chrom.setNativeID("");
      all_chroms.push_back(chrom);
    }
  }

  // cache parsed mz pairs
  std::vector<std::pair<double,double>> chrom_mz_cache(all_chroms.size(), std::make_pair(-1.0, -1.0));
  for (Size i = 0; i < all_chroms.size(); ++i)
  {
    String orig = all_chroms[i].getMetaValue("original_native_id");
    double c_prec = -1.0, c_prod = -1.0;
    if (!orig.empty())
    {
      if (sscanf(orig.c_str(), "%*[^0-9+\\-]%lf%*[^0-9+\\-]%lf", &c_prec, &c_prod) >= 2)
      {
        chrom_mz_cache[i] = std::make_pair(c_prec, c_prod);
      }
    }
  }

  // determine tolerances from cp
  double abs_tol = 0.5;
  double ppm_tol = 20.0;
  try
  {
    double cfg = cp.mz_extraction_window;
    String unit = cp.ppm ? "ppm" : "Th";
    if (cfg > 0)
    {
      if (unit == "ppm") { ppm_tol = cfg / 2.0; abs_tol = -1.0; }
      else { abs_tol = cfg / 2.0; ppm_tol = -1.0; }
    }
  }
  catch (...) {}

  // Map first matching transition to chromatogram
  for (const auto& tr : transition_exp.getTransitions())
  {
    double t_prec = tr.getPrecursorMZ();
    double t_prod = tr.getProductMZ();
    for (Size j = 0; j < chrom_mz_cache.size(); ++j)
    {
      if (!all_chroms[j].getNativeID().empty()) continue; // already assigned
      double c_prec = chrom_mz_cache[j].first;
      double c_prod = chrom_mz_cache[j].second;
      if (c_prec < 0 || c_prod < 0) continue;
      bool prec_ok = false;
      bool prod_ok = false;
      if (abs_tol >= 0) { prec_ok = (fabs(c_prec - t_prec) <= abs_tol); }
      if (!prec_ok && ppm_tol >= 0) { prec_ok = (fabs(c_prec - t_prec) / t_prec * 1e6 <= ppm_tol); }
      if (abs_tol >= 0) { prod_ok = (fabs(c_prod - t_prod) <= abs_tol); }
      if (!prod_ok && ppm_tol >= 0) { prod_ok = (fabs(c_prod - t_prod) / t_prod * 1e6 <= ppm_tol); }
      if (prec_ok && prod_ok)
      {
        all_chroms[j].setNativeID(tr.getNativeID());
        all_chroms[j].setMetaValue("mapped_transition", tr.getNativeID());
      }
    }
  }

  // filter and return mapped chromatograms
  std::vector<MSChromatogram> filtered;
  filtered.reserve(all_chroms.size());
  for (auto & c : all_chroms)
  {
    if (!c.getNativeID().empty()) filtered.push_back(c);
  }
  OPENMS_LOG_INFO << "SRM: mapped and will process " << filtered.size() << " chromatograms (from " << all_chroms.size() << ")" << std::endl;
  return filtered;
}

} // namespace OpenMS
