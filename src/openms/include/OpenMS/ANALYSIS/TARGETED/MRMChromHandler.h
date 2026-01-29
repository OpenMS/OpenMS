// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
// Need ChromExtractParams (defined in OpenSwathWorkflow.h)
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

// LightTargetedExperiment is defined in OpenSwath (openswathalgo). Forward-declare here to avoid
// including a non-existent OpenMS header location.
namespace OpenSwath { struct LightTargetedExperiment; }


namespace OpenMS
{
/**
  @brief Default SRM/MRM chromatogram provider declaration

  The `MRMChromHandler` class is the default implementation of
  `IChromatogramHandler` used by OpenSwathWorkflow when running in SRM/MRM
  mode. It delegates to the internal SRM/MRM helpers implemented in
  `MRMChromHandler.cpp` for chromatogram collection and mapping.

  This header exposes the provider's public interface (only the class
  declaration). The implementation is colocated in `MRMChromHandler.cpp` so the
  provider can remain a lightweight adapter while keeping implementation
  details private to the translation unit.
*/
class OPENMS_DLLAPI MRMChromHandler : public IChromatogramHandler
{
public:
  MRMChromHandler();
  ~MRMChromHandler() override;

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

  /// Utility function to normalize chromatogram precursor/product m/z and nativeID
  /// This handles chromatogram-only files where metadata may not be properly set
  static void normalizeChromatogramMZ(MSChromatogram& chrom);

};

} // namespace OpenMS
