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

  /**
    @brief Convenience entry point that runs @ref MRMFeatureFinderScoring on
           a list of already-extracted chromatograms.

    Intended for callers that already hold a list of MS chromatograms (XICs
    extracted upstream) and only need the feature-finding / scoring step
    without setting up a full SWATH-map / RT-transformation stack. Typical
    use is inside RT-normalisation flows in the OpenSwath workflow, where
    iRT chromatograms are extracted before the analytical-transition scoring
    runs.

    All scoring semantics come from @ref MRMFeatureFinderScoring; see its
    documentation for the parameter set passed in via
    @p feature_finder_param.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI ChromatogramProcessor
  {
  public:
    /// Default constructor (the class is stateless; the API is a single static method)
    ChromatogramProcessor() = default;

    /// Destructor
    ~ChromatogramProcessor() = default;

    /**
      @brief Score the supplied chromatograms against the transition list
             and emit features + transition groups.

      Tolerates transitions that have no matching chromatogram (they are
      simply skipped rather than reported as an error). Scoring is performed
      with no SWATH map and no RT transformation, so this overload is intended
      for callers that have already extracted the chromatograms upstream and
      only need the feature-finding step.

      Both outputs are populated in place: @p featureFile gets one feature
      per scored precursor, and @p transition_group_map gets one entry per
      transition group.

      @param[in]  chromatograms         Already-extracted XICs to score (one per transition).
      @param[in]  transition_exp        Transition list to align the chromatograms to.
      @param[in]  feature_finder_param  Parameter set forwarded to @ref MRMFeatureFinderScoring.
      @param[out] featureFile           Receives the scored features (one per precursor); not cleared first.
      @param[out] transition_group_map  Receives the per-transition-group scoring state.
    */
    static void pickExperiment(
      const std::vector<MSChromatogram> & chromatograms,
      const OpenSwath::LightTargetedExperiment & transition_exp,
      const Param & feature_finder_param,
      FeatureMap & featureFile,
      OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map);
  };

} // namespace OpenMS
