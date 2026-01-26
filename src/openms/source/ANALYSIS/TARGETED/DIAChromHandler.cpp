// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// 
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/DIAChromHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>

namespace OpenMS
{

std::vector<MSChromatogram> DIAChromHandler::collectIrtChromatogramsForIrt(
  const std::vector< OpenSwath::SwathMap > & swath_maps,
  const OpenSwath::LightTargetedExperiment & irt_transitions,
  const Param & mrm_mapping_param,
  const ChromExtractParams & cp,
  bool /*pasef*/,
  bool load_into_memory)
{
  std::vector<MSChromatogram> collected;

  ChromatogramExtractor extractor;
  for (Size i = 0; i < swath_maps.size(); ++i)
  {
    if (swath_maps[i].ms1) continue; // skip MS1 maps

    OpenSwath::SpectrumAccessPtr current = swath_maps[i].sptr;
    if (!current) continue;
    OpenSwath::SpectrumAccessPtr use_map = current;
    if (load_into_memory)
    {
      use_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>(new SpectrumAccessOpenMSInMemory(*current));
    }

  // Prepare coordinates for extraction of the iRT transitions
  std::vector< OpenSwath::ChromatogramPtr > chrom_list;
  std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
  // ChromatogramExtractor expects a non-const experiment type for template specialization,
  // so make a local copy to avoid instantiating the const specialization (which can cause
  // link errors for certain inline template specializations).
  OpenSwath::LightTargetedExperiment transition_exp_used;
  OpenSwathHelper::selectSwathTransitions(irt_transitions, transition_exp_used,
    cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
  ChromatogramExtractor::prepare_coordinates(chrom_list, coordinates, transition_exp_used, cp.rt_extraction_window, /*ms1*/ false, /*ms1_isotopes*/ 0);

    try
    {
      extractor.extractChromatograms(use_map, chrom_list, coordinates, cp.mz_extraction_window, cp.ppm, cp.im_extraction_window, cp.extraction_function);
      PeakMap pm;
      extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used, SpectrumSettings(), pm.getChromatograms(), false, cp.im_extraction_window);

      // append to collected
      for (auto & c : pm.getChromatograms())
      {
        collected.push_back(c);
      }
      OPENMS_LOG_DEBUG << "DIAChromHandler: extracted " << pm.getChromatograms().size() << " chromatograms from swath " << i << std::endl;
    }
    catch (const std::exception & e)
    {
      OPENMS_LOG_DEBUG << "DIAChromHandler: extraction failed for swath " << i << ": " << e.what() << std::endl;
      continue;
    }
  }

  // The ChromatogramExtractor::return_chromatogram call above already sets
  // the native ID and meta information for the extracted chromatograms
  // (including the transition native IDs). Remapping using MRMMapping is
  // redundant for DIA/SWATH-extracted chromatograms and can mis-assign
  // chromatograms when fallbacks are used. Return the extracted chromatograms
  // directly.
  OPENMS_LOG_DEBUG << "DIAChromHandler: returning " << collected.size() << " extracted iRT chromatograms" << std::endl;
  return collected;
}


std::vector<MSChromatogram> DIAChromHandler::extractAndMapChromatogramsForTransitions(
  const std::vector< OpenSwath::SwathMap > & swath_maps,
  const OpenSwath::LightTargetedExperiment & transition_exp,
  const ChromExtractParams & cp,
  const Param & mrm_mapping_param)
{
  std::vector<MSChromatogram> collected;
  ChromatogramExtractor extractor;

  for (Size i = 0; i < swath_maps.size(); ++i)
  {
    if (swath_maps[i].ms1) continue;
    OpenSwath::SpectrumAccessPtr current = swath_maps[i].sptr;
    if (!current) continue;

  std::vector< OpenSwath::ChromatogramPtr > chrom_list;
  std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
  OpenSwath::LightTargetedExperiment transition_exp_copy = transition_exp;
  ChromatogramExtractor::prepare_coordinates(chrom_list, coordinates, transition_exp_copy, cp.rt_extraction_window, /*ms1*/ false, /*ms1_isotopes*/ 0);

    try
    {
      extractor.extractChromatograms(current, chrom_list, coordinates, cp.mz_extraction_window, cp.ppm, cp.im_extraction_window, cp.extraction_function);
      PeakMap pm;
      extractor.return_chromatogram(chrom_list, coordinates, transition_exp, SpectrumSettings(), pm.getChromatograms(), false, cp.im_extraction_window);

  for (auto & c : pm.getChromatograms()) collected.push_back(c);
  OPENMS_LOG_DEBUG << "DIAChromHandler: extracted " << pm.getChromatograms().size() << " chromatograms from swath " << i << " for transitions" << std::endl;
    }
    catch (const std::exception & e)
    {
      OPENMS_LOG_DEBUG << "DIAChromHandler: extraction failed for swath " << i << " : " << e.what() << std::endl;
      continue;
    }
  }

  // ChromatogramExtractor::return_chromatogram already populated the
  // extracted chromatograms with native IDs corresponding to the transition
  // native IDs (coord.id). Therefore, we can directly return the collected
  // chromatograms without additional mapping.
  OPENMS_LOG_DEBUG << "DIAChromHandler: returning " << collected.size() << " extracted chromatograms for transitions" << std::endl;
  return collected;
}

} // namespace OpenMS
