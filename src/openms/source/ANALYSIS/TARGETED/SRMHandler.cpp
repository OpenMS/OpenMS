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
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>

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
      // copy precursor/product information from the underlying OpenMS chromatogram metadata
      ChromatogramSettings cs = openms_map->getChromatogramMetaInfo(i);
      chrom.setPrecursor(cs.getPrecursor());
      chrom.setProduct(cs.getProduct());
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
  const ChromExtractParams & cp,
  const Param & mrm_mapping_param)
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
      // copy precursor/product information from the underlying OpenMS chromatogram metadata
      ChromatogramSettings cs = openms_map->getChromatogramMetaInfo(i);
      chrom.setPrecursor(cs.getPrecursor());
      chrom.setProduct(cs.getProduct());
      String orig = openms_map->getChromatogramNativeID(i);
      chrom.setMetaValue("original_native_id", orig);
      // clear nativeID until mapped
      chrom.setNativeID("");
      all_chroms.push_back(chrom);
    }
  }

  // Use canonical MRMMapping to map chromatograms to transitions by default.
  // Build a PeakMap from the collected MSChromatogram objects and call MRMMapping.
  // Quick diagnostics: count how many chromatograms have precursor/product m/z present
  Size prec_set_count = 0, prod_set_count = 0, both_set_count = 0;
  for (const auto & c : all_chroms)
  {
    bool pset = std::fabs(c.getPrecursor().getMZ()) > 1e-8;
    bool oset = std::fabs(c.getProduct().getMZ()) > 1e-8;
    if (pset) ++prec_set_count;
    if (oset) ++prod_set_count;
    if (pset && oset) ++both_set_count;
  }
  OPENMS_LOG_INFO << "SRM: chromatograms with precursor m/z=" << prec_set_count << ", product m/z=" << prod_set_count << ", both=" << both_set_count << " (total=" << all_chroms.size() << ")" << std::endl;

  PeakMap input_pm;
  input_pm.getChromatograms() = all_chroms; // copy all chromatograms (meta-values like original_native_id preserved)

  PeakMap output_pm;
  MRMMapping mapper;
  // Apply user-supplied mapping parameters (e.g., map_multiple_assays)
  try
  {
    mapper.setParameters(mrm_mapping_param);
    mapper.mapExperiment(input_pm, transition_exp, output_pm);
  }
  catch (const std::exception & e)
  {
    OPENMS_LOG_WARN << "SRM: MRMMapping threw an exception: " << e.what() << ". Falling back to no mappings." << std::endl;
    std::vector<MSChromatogram> empty;
    return empty;
  }

  // Collect mapped chromatograms (those that received a nativeID by the mapper)
  std::vector<MSChromatogram> filtered;
  filtered.reserve(output_pm.getChromatograms().size());
  Size mapped_count = 0;
  for (const auto & c : output_pm.getChromatograms())
  {
    if (!c.getNativeID().empty())
    {
      ++mapped_count;
      filtered.push_back(c);
    }
  }
  Size total = all_chroms.size();
  OPENMS_LOG_INFO << "SRM: mapped and will process " << mapped_count << " chromatograms (from " << total << ")" << std::endl;
  if (mapped_count > 0)
  {
    return filtered;
  }

  // Fallback: if MRMMapping produced no mappings, try isolation-window based mapping
  OPENMS_LOG_INFO << "SRM: MRMMapping produced 0 mappings, falling back to isolation-window based mapping." << std::endl;

  // cache parsed mz pairs AND isolation-window info from chromatograms
  // For each chromatogram we record:
  //  - parsed numeric precursor/product (fallback when isolation windows are not present)
  //  - precursor/product isolation window ranges (when present in the chromatogram metadata)
  struct ChromInfo { double prec_lower{0.0}; double prec_upper{0.0}; bool has_prec_win{false};
                     double prod_lower{0.0}; double prod_upper{0.0}; bool has_prod_win{false}; };

  std::vector<ChromInfo> chrom_info(all_chroms.size());
  for (Size i = 0; i < all_chroms.size(); ++i)
  {
    // Use explicit precursor/product isolation window information from the chromatogram settings only.
    const Precursor &prec = all_chroms[i].getPrecursor();
    const Product &prod = all_chroms[i].getProduct();
    if (prec.getMZ() > 0.0)
    {
      chrom_info[i].prec_lower = prec.getMZ() - prec.getIsolationWindowLowerOffset();
      chrom_info[i].prec_upper = prec.getMZ() + prec.getIsolationWindowUpperOffset();
      chrom_info[i].has_prec_win = true;
    }
    if (prod.getMZ() > 0.0)
    {
      chrom_info[i].prod_lower = prod.getMZ() - prod.getIsolationWindowLowerOffset();
      chrom_info[i].prod_upper = prod.getMZ() + prod.getIsolationWindowUpperOffset();
      chrom_info[i].has_prod_win = true;
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

  // Map first matching transition to chromatogram preferring explicit isolation windows
  for (const auto& tr : transition_exp.getTransitions())
  {
    double t_prec = tr.getPrecursorMZ();
    double t_prod = tr.getProductMZ();
    for (Size j = 0; j < chrom_info.size(); ++j)
    {
      if (!all_chroms[j].getNativeID().empty()) continue; // already assigned

      bool prec_ok = false;
      bool prod_ok = false;

      // 1) Try isolation-window based matching when present in the chromatogram
      if (chrom_info[j].has_prec_win)
      {
        double lower = chrom_info[j].prec_lower;
        double upper = chrom_info[j].prec_upper;
        if (abs_tol >= 0)
        {
          if ((t_prec >= lower - abs_tol) && (t_prec <= upper + abs_tol)) prec_ok = true;
        }
        if (!prec_ok && ppm_tol >= 0)
        {
          double ext = t_prec * ppm_tol * 1e-6;
          if ((t_prec >= lower - ext) && (t_prec <= upper + ext)) prec_ok = true;
        }
      }

      if (chrom_info[j].has_prod_win)
      {
        double lower = chrom_info[j].prod_lower;
        double upper = chrom_info[j].prod_upper;
        if (abs_tol >= 0)
        {
          if ((t_prod >= lower - abs_tol) && (t_prod <= upper + abs_tol)) prod_ok = true;
        }
        if (!prod_ok && ppm_tol >= 0)
        {
          double ext = t_prod * ppm_tol * 1e-6;
          if ((t_prod >= lower - ext) && (t_prod <= upper + ext)) prod_ok = true;
        }
      }

      // Only map when BOTH precursor and product isolation windows are present and matched.
      if (chrom_info[j].has_prec_win && chrom_info[j].has_prod_win && prec_ok && prod_ok)
      {
        all_chroms[j].setNativeID(tr.getNativeID());
        all_chroms[j].setMetaValue("mapped_transition", tr.getNativeID());
        filtered.push_back(all_chroms[j]);
      }
    }
  }

  OPENMS_LOG_INFO << "SRM (fallback): mapped and will process " << filtered.size() << " chromatograms (from " << all_chroms.size() << ")" << std::endl;
  return filtered;
}

} // namespace OpenMS
