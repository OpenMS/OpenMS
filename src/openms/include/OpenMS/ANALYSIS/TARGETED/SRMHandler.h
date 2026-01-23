// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
// Need ChromExtractParams (defined in OpenSwathWorkflow.h)
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>

// LightTargetedExperiment is defined in OpenSwath (openswathalgo). Forward-declare here to avoid
// including a non-existent OpenMS header location.
namespace OpenSwath { struct LightTargetedExperiment; }

namespace OpenMS
{
class SRMHandler
{
public:
  /// Collect chromatograms from chromatogram-only swath maps and return them with original nativeIDs preserved
  static std::vector<OpenMS::MSChromatogram> collectChromatogramsForIrt(const std::vector< OpenSwath::SwathMap > & swath_maps);

  /// From chromatogram-only swath maps, map chromatograms to transitions (by nativeID exact match or mz-based fallback)
  /// and return only the chromatograms that were successfully mapped (their nativeID will be set to the matched transition nativeID).
  static std::vector<OpenMS::MSChromatogram> extractAndMapChromatogramsForTransitions(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment & transition_exp,
    const ChromExtractParams & cp);
};

} // namespace OpenMS
