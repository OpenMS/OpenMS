// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/SRMHandler.h>
#include <OpenMS/ANALYSIS/TARGETED/DefaultChromHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h>
#include <OpenMS/KERNEL/MSChromatogram.h>


namespace OpenMS
{

  SRMChromHandler::SRMChromHandler() = default;
  SRMChromHandler::~SRMChromHandler() = default;

  std::vector<MSChromatogram> SRMChromHandler::collectIrtChromatogramsForIrt(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & irt_transitions,
    const Param & mrm_mapping_param,
    const ChromExtractParams & /*cp*/,
    bool /*pasef*/,
    bool /*load_into_memory*/)
  {
    std::vector<MSChromatogram> irt_chroms;
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
        irt_chroms.push_back(chrom);
        ++total_chroms;
      }
    }
    OPENMS_LOG_DEBUG << "SRMChromHandler: collected " << irt_chroms.size() << " raw iRT chromatograms." << std::endl;

    // Try mapping using MRMMapping with provided parameters
    try
    {
      PeakMap pm;
      pm.setChromatograms(irt_chroms);
      OPENMS_LOG_DEBUG << "SRMChromHandler: PeakMap prepared with " << pm.getChromatograms().size() << " chromatograms." << std::endl;
      PeakMap mapped_pm;
      MRMMapping mapper;
      Param mrm_p = mapper.getParameters();
      mrm_p.update(mrm_mapping_param, false, false, false, false, getGlobalLogDebug());
      mrm_p.setValue("map_multiple_assays", "true");
      mapper.setParameters(mrm_p);
      mapper.mapExperiment(pm, irt_transitions, mapped_pm);
      std::vector<MSChromatogram> result = mapped_pm.getChromatograms();
      OPENMS_LOG_DEBUG << "SRMChromHandler: mapping produced " << result.size() << " chromatograms (mapped)." << std::endl;
      return result;
    }
    catch (const std::exception & e)
    {
      // Mapping failed for iRT chromatograms: this is fatal for downstream
      // processing (can't calibrate without iRT matches). 
      OPENMS_LOG_ERROR << "SRMChromHandler: MRMMapping failed: " << e.what()
                       << " - aborting (iRT mapping required)." << std::endl;
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String("MRMMapping failed during iRT collection: ") + e.what());
    }
  }

  std::vector<MSChromatogram> SRMChromHandler::extractAndMapChromatogramsForTransitions(
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
    Size prec_set_count = 0, prod_set_count = 0, both_set_count = 0;
    for (const auto & c : all_chroms)
    {
      bool pset = std::fabs(c.getPrecursor().getMZ()) > 1e-8;
      bool oset = std::fabs(c.getProduct().getMZ()) > 1e-8;
      if (pset) ++prec_set_count;
      if (oset) ++prod_set_count;
      if (pset && oset) ++both_set_count;
    }
    OPENMS_LOG_DEBUG << "SRM: chromatograms with precursor m/z=" << prec_set_count << ", product m/z=" << prod_set_count << ", both=" << both_set_count << " (total=" << all_chroms.size() << ")" << std::endl;

    PeakMap input_pm;
    input_pm.getChromatograms() = all_chroms; 

    PeakMap output_pm;
    MRMMapping mapper;
    try
    {
      mapper.setParameters(mrm_mapping_param);
      mapper.mapExperiment(input_pm, transition_exp, output_pm);
    }
    catch (const std::exception & e)
    {
      OPENMS_LOG_ERROR << "SRMChromHandler: MRMMapping failed: " << e.what()
                      << " - aborting (transition-chromatogram experiment mapping)." << std::endl;
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        String("MRMMapping failed during transition-chromatogram experiment mapping collection: ") + e.what());
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
    OPENMS_LOG_DEBUG << "SRM: mapped and will process " << mapped_count << " chromatograms (from " << total << ")" << std::endl;
    if (mapped_count > 0)
    {
      OPENMS_LOG_DEBUG << "SRMChromHandler: SRMHandler returned " << filtered.size() << " chromatograms." << std::endl;
      return filtered;
    }
  }

  std::unique_ptr<IChromatogramHandler> IChromatogramHandler::createDefault()
  {
    // Return the default delegating handler which will pick SRM vs DIA at runtime
    return std::unique_ptr<IChromatogramHandler>(new DefaultChromHandler());
  }

} // namespace OpenMS
