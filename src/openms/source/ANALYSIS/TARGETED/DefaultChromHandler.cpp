// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// 
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/DefaultChromHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  DefaultChromHandler::DefaultChromHandler()
    : srm_(new SRMChromHandler()), dia_(new DIAChromHandler())
  {
  }

  DefaultChromHandler::~DefaultChromHandler() = default;

  std::vector<MSChromatogram> DefaultChromHandler::collectIrtChromatogramsForIrt(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & irt_transitions,
    const Param & mrm_mapping_param,
    const ChromExtractParams & cp,
    bool pasef,
    bool load_into_memory)
  {
    // Decide whether inputs are chromatogram-only (SRM) or spectral (DIA)
    bool srm_mode = true;
    for (const auto & sm : swath_maps)
    {
      if (sm.ms1 || sm.sptr->getNrSpectra() > 0)
      {
        srm_mode = false;
        break;
      }
    }

    if (srm_mode)
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating iRT collection to SRM handler" << std::endl;
      return srm_->collectIrtChromatogramsForIrt(swath_maps, irt_transitions, mrm_mapping_param, cp, pasef, load_into_memory);
    }
    else
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating iRT collection to DIA handler" << std::endl;
      return dia_->collectIrtChromatogramsForIrt(swath_maps, irt_transitions, mrm_mapping_param, cp, pasef, load_into_memory);
    }
  }

  std::vector<MSChromatogram> DefaultChromHandler::extractAndMapChromatogramsForTransitions(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & transition_exp,
    const ChromExtractParams & cp,
    const Param & mrm_mapping_param)
  {
    bool srm_mode = true;
    for (const auto & sm : swath_maps)
    {
      if (sm.ms1 || sm.sptr->getNrSpectra() > 0)
      {
        srm_mode = false;
        break;
      }
    }

    if (srm_mode)
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating transition extraction to SRM handler" << std::endl;
      return srm_->extractAndMapChromatogramsForTransitions(swath_maps, transition_exp, cp, mrm_mapping_param);
    }
    else
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating transition extraction to DIA handler" << std::endl;
      return dia_->extractAndMapChromatogramsForTransitions(swath_maps, transition_exp, cp, mrm_mapping_param);
    }
  }

} // namespace OpenMS
