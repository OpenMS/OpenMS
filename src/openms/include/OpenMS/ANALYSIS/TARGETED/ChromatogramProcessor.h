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
    @brief Thin wrapper around MRMFeatureFinderScoring::pickExperiment for already-extracted chromatograms.

    The class exists so callers that already hold a list of MS chromatograms (i.e. XICs that
    have been extracted upstream) can drive feature picking / scoring without having to
    construct an empty SWATH-map / transformation stack themselves. It is typically used
    inside RT-normalisation flows in the OpenSwath workflow, where iRT chromatograms have
    been extracted before the analytical-transition scoring runs.

    All work is delegated to @ref MRMFeatureFinderScoring; see its documentation for the
    semantics of the parameter set passed in through @p feature_finder_param.

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
      @brief Score the supplied chromatograms against the transition list and emit features + transition groups.

      Wraps @p chromatograms in a @ref PeakMap and a @ref SpectrumAccessOpenMS accessor,
      constructs an @ref MRMFeatureFinderScoring instance parameterised with
      @p feature_finder_param (and with @c setStrictFlag(false) so the call does not abort
      on transitions that have no matching chromatogram), then calls @c pickExperiment()
      with an empty SWATH-map vector and an identity @ref TransformationDescription.

      Both outputs are populated in place: @p featureFile gets one feature per scored
      precursor, and @p transition_group_map gets one entry per transition group.

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
