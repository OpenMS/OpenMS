// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/KERNEL/MSChromatogram.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <memory>

// Forward-declare OpenSwath types at global namespace (defined in openswathalgo)
namespace OpenSwath { struct LightTargetedExperiment; struct SwathMap; }

namespace OpenMS
{
  // ChromExtractParams is defined in OpenSwathWorkflow.h (OpenMS namespace)
  struct ChromExtractParams;

  /**
    @brief Abstract interface for providing chromatograms 

    This interface isolates the rest of the OpenSwath workflow from
    the source and mapping strategy for chromatograms. Implementations are
    expected to provide two main services:

    - collectIrtChromatogramsForIrt(...): collect candidate chromatograms used
      for iRT calibration and attempt to map them to the supplied iRT
      transitions. 

    - extractAndMapChromatogramsForTransitions(...): either select existing
      SRM/MRM chromatograms or extract chromatograms from full-scan data
      (DIA/SWATH/PRM) and map them to the supplied transitions. Implementations
      should return only chromatograms that are mapped and ready for
      processing.
  */
  class OPENMS_DLLAPI IChromatogramHandler
  {
  public:
    IChromatogramHandler() = default;
    virtual ~IChromatogramHandler() = default;

    /// Collect iRT chromatograms from the swath maps and try to map them to the provided iRT transitions
    virtual std::vector<MSChromatogram> collectIrtChromatogramsForIrt(
      const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenSwath::LightTargetedExperiment & irt_transitions,
      const Param & mrm_mapping_param,
      const ChromExtractParams & cp,
      const TransformationDescription& trafo = TransformationDescription(),
      bool pasef = false,
      bool load_into_memory = false) = 0;

    /// Extract (or select) chromatograms for the given transitions and return mapped & filtered chromatograms
    virtual std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(
      const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenSwath::LightTargetedExperiment & transition_exp,
      const ChromExtractParams & cp,
      const Param & mrm_mapping_param) = 0;

    /// Factory: create the default handler (currently SRM/MRM-based)
    static std::unique_ptr<IChromatogramHandler> createDefault();
  };

} // namespace OpenMS
