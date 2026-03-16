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
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  DefaultChromHandler::DefaultChromHandler()
    : mrm_(new MRMChromHandler()), dia_(new DIAChromHandler())
  {
  }

  DefaultChromHandler::~DefaultChromHandler() = default;

  std::vector<MSChromatogram> DefaultChromHandler::collectIrtChromatogramsForIrt(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & irt_transitions,
    const Param & mrm_mapping_param,
    const ChromExtractParams & cp,
    const TransformationDescription& trafo,
    bool pasef,
    bool load_into_memory)
  {
    // Decide whether inputs are chromatogram-only (SRM/MRM) or spectral (DIA)
    bool mrm_mode = true;
    for (const auto & sm : swath_maps)
    {
      if (sm.ms1 || (sm.sptr && sm.sptr->getNrSpectra() > 0))
      {
        mrm_mode = false;
        break;
      }
    }

    if (mrm_mode)
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating iRT collection to SRM/MRM handler" << std::endl;
      try
      {
        return mrm_->collectIrtChromatogramsForIrt(swath_maps, irt_transitions, mrm_mapping_param, cp, trafo, pasef, load_into_memory);
      }
      catch (const Exception::IllegalArgument& e)
      {
        // No chromatograms mapped to transitions: treat as a recoverable condition for SRM/MRM data.
        // SRM/MRM mzML inputs are chromatogram-only (no spectra) and mapping relies on mapping precursor/product m/z that may not match the provided transition list. In these cases there simply are no chroms to return.
        // We remove the empty chroms downstream and only process those that were mapped successfully.
        OPENMS_LOG_WARN << "DefaultChromHandler: SRM/MRM handler reported no iRT chromatograms: " << e.what() << " - returning empty result" << std::endl;
        return std::vector<MSChromatogram>();
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_ERROR << "DefaultChromHandler: SRM/MRM handler failed for iRT collection: " << e.what() << std::endl;
        throw;
      }
    }
    else
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating iRT collection to DIA handler" << std::endl;
      try
      {
        return dia_->collectIrtChromatogramsForIrt(swath_maps, irt_transitions, mrm_mapping_param, cp, trafo, pasef, load_into_memory);
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_ERROR << "DefaultChromHandler: DIA handler failed for iRT collection: " << e.what() << std::endl;
        throw;
      }
    }
  }

  std::vector<MSChromatogram> DefaultChromHandler::extractAndMapChromatogramsForTransitions(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & transition_exp,
    const ChromExtractParams & cp,
    const Param & mrm_mapping_param)
  {
    bool mrm_mode = true;
    for (const auto & sm : swath_maps)
    {
      if (sm.ms1 || sm.sptr->getNrSpectra() > 0)
      {
        mrm_mode = false;
        break;
      }
    }

    if (mrm_mode)
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating transition extraction to SRM/MRM handler" << std::endl;
      try
      {
        return mrm_->extractAndMapChromatogramsForTransitions(swath_maps, transition_exp, cp, mrm_mapping_param);
      }
      catch (const Exception::IllegalArgument& e)
      {
        // No chromatograms mapped to transitions: treat as a recoverable condition for SRM/MRM data.
        // SRM/MRM mzML inputs are chromatogram-only (no spectra) and mapping relies on mapping precursor/product m/z that may not match the provided transition list. In these cases there simply are no chroms to return.
        // We remove the empty chroms downstream and only process those that were mapped successfully.
        OPENMS_LOG_WARN << "DefaultChromHandler: SRM/MRM handler reported no chromatograms mapped to transitions: " << e.what() << " - returning empty result" << std::endl;
        return std::vector<MSChromatogram>();
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_ERROR << "DefaultChromHandler: SRM/MRM handler failed for transition extraction: " << e.what() << std::endl;
        throw;
      }
    }
    else
    {
      OPENMS_LOG_DEBUG << "DefaultChromHandler: delegating transition extraction to DIA handler" << std::endl;
      try
      {
        return dia_->extractAndMapChromatogramsForTransitions(swath_maps, transition_exp, cp, mrm_mapping_param);
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_ERROR << "DefaultChromHandler: DIA handler failed for transition extraction: " << e.what() << std::endl;
        throw;
      }
    }
  }

} // namespace OpenMS
