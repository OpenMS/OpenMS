// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ML/RANSAC/RANSACModel.h>

namespace OpenMS
{

  namespace Math
  {
    /**
      @brief Implementation of a quadratic RANSAC model fit.
      
      Using generic plug-in template base class 'RansacModel' using 'Curiously recurring template pattern' (CRTP).
    */
    class OPENMS_DLLAPI RansacModelQuadratic
      : public RansacModel<RansacModelQuadratic>
    {
    public:
      static ModelParameters rm_fit_impl(const DVecIt& begin, const DVecIt& end);

      static double rm_rsq_impl(const DVecIt& begin, const DVecIt& end);

      static double rm_rss_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients);

      static DVec rm_inliers_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold);

    };

  }


}
