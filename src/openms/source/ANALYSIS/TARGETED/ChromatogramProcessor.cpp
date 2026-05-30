// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/ChromatogramProcessor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>

namespace OpenMS
{

void ChromatogramProcessor::pickExperiment(
  const std::vector<MSChromatogram> & chromatograms,
  const OpenSwath::LightTargetedExperiment & transition_exp,
  const Param & feature_finder_param,
  FeatureMap & featureFile,
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map)
{
  OPENMS_LOG_INFO << "ChromatogramProcessor::pickExperiment() called with " << chromatograms.size() << " chromatograms." << std::endl;
  MRMFeatureFinderScoring featureFinder;
  const Param& ffparam(feature_finder_param);
  featureFinder.setParameters(ffparam);

  // Prepare the data with the chromatograms
  std::shared_ptr<PeakMap > xic_map(new PeakMap);
  xic_map->setChromatograms(chromatograms);
  OPENMS_LOG_INFO << "ChromatogramProcessor: PeakMap populated with " << xic_map->getChromatograms().size() << " chromatograms." << std::endl;
  OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(xic_map));

  // Call pickExperiment which will fill featureFile and transition_group_map
  std::vector<OpenSwath::SwathMap> empty_swath_maps;
  TransformationDescription empty_trafo;
  featureFinder.setStrictFlag(false);
  OPENMS_LOG_INFO << "ChromatogramProcessor: invoking MRMFeatureFinderScoring::pickExperiment()" << std::endl;
  featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp, empty_trafo, empty_swath_maps, transition_group_map);
  OPENMS_LOG_INFO << "ChromatogramProcessor: pickExperiment() completed. Extracted " << transition_group_map.size() << " transition groups." << std::endl;
}

} // namespace OpenMS
