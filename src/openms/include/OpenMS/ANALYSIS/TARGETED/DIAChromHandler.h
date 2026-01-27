// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// 
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>

namespace OpenMS
{
  class OPENMS_DLLAPI DIAChromHandler : public IChromatogramHandler
  {
  public:
    DIAChromHandler() = default;
    ~DIAChromHandler() override = default;

    std::vector<MSChromatogram> collectIrtChromatogramsForIrt(
      const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenSwath::LightTargetedExperiment & irt_transitions,
      const Param & mrm_mapping_param,
      const ChromExtractParams & cp,
      const TransformationDescription& trafo = TransformationDescription(),
      bool pasef = false,
      bool load_into_memory = false) override;

    std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(
      const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenSwath::LightTargetedExperiment & transition_exp,
      const ChromExtractParams & cp,
      const Param & mrm_mapping_param) override;
  };

} // namespace OpenMS
