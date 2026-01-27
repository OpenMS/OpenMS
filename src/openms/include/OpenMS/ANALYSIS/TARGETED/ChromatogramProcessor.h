// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/KERNEL/MSChromatogram.h>
// forward-declared below
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>

// Forward-declare OpenSwath types at global namespace (defined in openswathalgo)
namespace OpenSwath { struct LightTargetedExperiment; }

namespace OpenMS
{

  /// Small helper to encapsulate feature finding / scoring on chromatogram inputs.
  /// This isolates MRMFeatureFinderScoring setup used in RT normalization and elsewhere.
  class OPENMS_DLLAPI ChromatogramProcessor
  {
  public:
    ChromatogramProcessor() = default;
    ~ChromatogramProcessor() = default;

    /// Run pickExperiment on the provided chromatograms and fill featureFile and transition_group_map
    static void pickExperiment(
      const std::vector<MSChromatogram> & chromatograms,
      const OpenSwath::LightTargetedExperiment & transition_exp,
      const Param & feature_finder_param,
      FeatureMap & featureFile,
      OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map);
  };

} // namespace OpenMS
