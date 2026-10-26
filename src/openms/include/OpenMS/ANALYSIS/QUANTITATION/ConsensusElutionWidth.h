// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

#include <vector>

namespace OpenMS
{
  class ConsensusMap;
  class FeatureMap;
  class TransformationDescription;

  /**
    @brief Re-quantifies linked features at one elution width per precursor

    FeatureFinderIdentification fits an elution model to each feature and reports its area, height x sigma x sqrt(2 pi)
    for the Gaussian model. The width is fitted per run, so its error moves the area of a precursor from run to run
    although the peak height is right. This class replaces each run's width by one consensus width per consensus
    feature:

      w_c = median over the identified features of the consensus feature of sigma_r * g_r'(t_r)

    where t_r is the feature's apex in run r (meta value "model_center", run RT) and g_r' the slope of the RT
    transformation of run r (run RT to consensus RT) at t_r, which carries the width onto the consensus RT scale.
    Every feature of the consensus feature with a valid fit then reports

      height_r x (w_c / g_r'(t_r)) x sqrt(2 pi).

    A feature is identified if it carries a peptide identification of its own run. Features without a valid Gaussian
    fit (meta values "model_height", "model_Gauss_sigma", "model_status" starting with '0') keep their intensity, and so
    do all features of a consensus feature without identified fitted features.

    Ported from bigbio/andes #150 (finalize: one consensus width per precursor).
  */
  class OPENMS_DLLAPI ConsensusElutionWidth
  {
  public:
    /// feature meta values of the elution model fit that apply() reads (written by ElutionModelFitter)
    static constexpr const char* FIT_META_VALUES[] = {"model_height", "model_Gauss_sigma", "model_center", "model_status"};

    struct OPENMS_DLLAPI Stats
    {
      Size consensus_features = 0;   ///< consensus features with a consensus width
      Size handles_updated = 0;      ///< feature handles re-quantified
      Size handles_kept = 0;         ///< feature handles that kept their intensity
    };

    /**
      @brief Re-quantify the feature handles of @p consensus

      @param consensus The linked map; its handle intensities are replaced
      @param feature_maps The feature maps that were linked, indexed by map index (features are looked up by unique id)
      @param transformations RT transformations (run RT to consensus RT) per map; empty or shorter: slope 1 for the rest
    */
    static Stats apply(ConsensusMap& consensus,
                       const std::vector<FeatureMap>& feature_maps,
                       const std::vector<TransformationDescription>& transformations);

    /// Slope of @p trafo at @p rt by a central difference over +- 1 s (1 for an identity or empty transformation)
    static double localSlope(const TransformationDescription& trafo, double rt);
  };
}
