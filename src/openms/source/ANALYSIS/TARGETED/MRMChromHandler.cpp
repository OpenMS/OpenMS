// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/MRMChromHandler.h>
#include <OpenMS/ANALYSIS/TARGETED/DefaultChromHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

#include <cstdio>


namespace OpenMS
{

  MRMChromHandler::MRMChromHandler() = default;
  MRMChromHandler::~MRMChromHandler() = default;

  void MRMChromHandler::normalizeChromatogramMZ(MSChromatogram& chrom)
  {
    double prec_mz = chrom.getPrecursor().getMZ();
    double prod_mz = chrom.getProduct().getMZ();

    // If precursor/product m/z are not set, try to read them from CV terms
    // (e.g. isolation window target m/z CV term MS:1000827) which some
    // mzML writers put into the <precursor>/<product> cvParams rather than
    // the Precursor/Product mz fields.
    if (!(prec_mz > 0))
    {
      const auto & prec_terms = chrom.getPrecursor().getCVTerms();
      auto it = prec_terms.find("MS:1000827");
      if (it != prec_terms.end() && !it->second.empty())
      {
        try
        {
          prec_mz = it->second[0].getValue().toString().toDouble();
          Precursor p = chrom.getPrecursor();
          p.setMZ(prec_mz);
          chrom.setPrecursor(p);
        }
        catch (const std::exception & e)
        {
          OPENMS_LOG_WARN << "Failed to parse precursor m/z for chromatogram (nativeID="
                          << chrom.getNativeID() << "): " << e.what() << std::endl;
        }
        catch (...)
        {
          OPENMS_LOG_WARN << "Unknown error while parsing precursor m/z for chromatogram (nativeID="
                          << chrom.getNativeID() << ")." << std::endl;
        }
      }
    }

    if (!(prod_mz > 0))
    {
      const auto & prod_terms = chrom.getProduct().getCVTerms();
      auto it = prod_terms.find("MS:1000827");
      if (it != prod_terms.end() && !it->second.empty())
      {
        prod_mz = it->second[0].getValue().toString().toDouble();
        Product q = chrom.getProduct();
        q.setMZ(prod_mz);
        chrom.setProduct(q);
      }
    }

    // Ensure meta-values and a canonical nativeID exist for whichever values we have
    if (prec_mz > 0)
    {
      chrom.setMetaValue("precursor_mz", prec_mz);
    }
    if (prod_mz > 0)
    {
      chrom.setMetaValue("product_mz", prod_mz);
    }

    // Build a compact canonical nativeID depending on available values
    String nid_out;
    if (prec_mz > 0 && prod_mz > 0)
    {
      nid_out = String("precursor=") + String(prec_mz) + ",product=" + String(prod_mz);
    }
    else if (prod_mz > 0)
    {
      nid_out = String("product=") + String(prod_mz);
    }
    else if (prec_mz > 0)
    {
      nid_out = String("precursor=") + String(prec_mz);
    }

    if (!nid_out.empty()) chrom.setNativeID(nid_out);
  }

  std::vector<MSChromatogram> MRMChromHandler::collectIrtChromatogramsForIrt(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & irt_transitions,
    const Param & mrm_mapping_param,
    const ChromExtractParams & /*cp*/,
    const TransformationDescription& /*trafo*/,
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
        
        // Normalize chromatogram m/z values using the utility function
        normalizeChromatogramMZ(chrom);
        
        // For iRT calibration, preserve the original nativeID in the chromatogram
        chrom.setNativeID(native_id);
        irt_chroms.push_back(chrom);
        ++total_chroms;
      }
    }

    // Try mapping using MRMMapping with provided parameters
    try
    {
      PeakMap pm;
      pm.setChromatograms(irt_chroms);
      OPENMS_LOG_DEBUG << "PeakMap prepared with " << pm.getChromatograms().size() << " chromatograms for MRMMapping." << std::endl;
      PeakMap mapped_pm;
      MRMMapping mapper;
      Param mrm_p = mapper.getParameters();
      mrm_p.update(mrm_mapping_param, false, false, false, false, getGlobalLogDebug());
      // For iRT calibration we prefer to allow multiple assay mappings to
      // maximize the number of candidate iRT chromatograms mapped. However,
      // preserve explicit user intent: if the user provided a value for
      // `map_multiple_assays` in `mrm_mapping_param`, respect that setting
      // instead of forcing an override.
      if (!mrm_mapping_param.exists("map_multiple_assays"))
      {
        mrm_p.setValue("map_multiple_assays", "true");
      }
      else
      {
        // Respect user-specified value; log at debug level for traceability.
        OPENMS_LOG_DEBUG << "MRMChromHandler: preserving user-provided 'map_multiple_assays' = "
                         << mrm_mapping_param.getValue("map_multiple_assays").toString()
                         << " for iRT mapping." << std::endl;
      }
      mapper.setParameters(mrm_p);
      mapper.mapExperiment(pm, irt_transitions, mapped_pm);
      std::vector<MSChromatogram> result = mapped_pm.getChromatograms();
      OPENMS_LOG_DEBUG << "MRMMapping produced " << result.size() << " chromatograms." << std::endl;
      
      // Filter out empty chromatograms
      std::vector<MSChromatogram> filtered_result;
      for (const auto& chrom : result) {
        if (!chrom.empty()) {
          filtered_result.push_back(chrom);
        }
      }
      OPENMS_LOG_INFO << "Mapped " << filtered_result.size() << " iRT chromatograms out of " << total_chroms << " total chromatograms." << std::endl;
      return filtered_result;
    }
    catch (const std::exception & e)
    {
      // Mapping failed for iRT chromatograms: this is fatal for downstream
      // processing (can't calibrate without iRT matches). 
      OPENMS_LOG_ERROR << "MRMChromHandler: MRMMapping failed: " << e.what()
                       << " - aborting (iRT mapping required)." << std::endl;
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String("MRMMapping failed during iRT collection: ") + e.what());
    }
  }

  std::vector<MSChromatogram> MRMChromHandler::extractAndMapChromatogramsForTransitions(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & transition_exp,
    const ChromExtractParams & cp,
    const Param & mrm_mapping_param)
  {
    // cp parameter is not used in SRM/MRM handler as it works with pre-existing chromatograms
    (void)cp;
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
        
        // Normalize chromatogram m/z values using the utility function
        normalizeChromatogramMZ(chrom);
        
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
    OPENMS_LOG_DEBUG << "SRM/MRM: chromatograms with precursor m/z=" << prec_set_count << ", product m/z=" << prod_set_count << ", both=" << both_set_count << " (total=" << all_chroms.size() << ")" << std::endl;

    // Filter out empty chromatograms before mapping
    std::vector<MSChromatogram> non_empty_chroms;
    for (const auto& chrom : all_chroms) {
      if (!chrom.empty()) {
        non_empty_chroms.push_back(chrom);
      }
    }

    PeakMap input_pm;
    input_pm.getChromatograms() = non_empty_chroms; 

    PeakMap output_pm;
    MRMMapping mapper;
    try
    {
      mapper.setParameters(mrm_mapping_param);
      mapper.mapExperiment(input_pm, transition_exp, output_pm);
    }
    catch (const std::exception & e)
    {
      OPENMS_LOG_ERROR << "MRMChromHandler: MRMMapping failed: " << e.what()
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
    if (mapped_count > 0)
    {
      OPENMS_LOG_INFO << "Mapped " << mapped_count << " chromatograms to transitions out of " << total << " total chromatograms." << std::endl;
      return filtered;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "There were no chromatograms mapped to the input transition list.");
    }
  }

  std::unique_ptr<IChromatogramHandler> IChromatogramHandler::createDefault()
  {
    // Return the default delegating handler which will pick SRM/MRM vs DIA at runtime
    return std::unique_ptr<IChromatogramHandler>(new DefaultChromHandler());
  }

} // namespace OpenMS
