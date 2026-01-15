// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  namespace Math
  {

    /**
      @brief Generic plug-in template base class using 'Curiously recurring template pattern' (CRTP)
             to allow for arbitrary RANSAC models (e.g. linear or quadratic fits).

      An actual model should derive from this class using CRTP (see RansacModelLinear for an example)
      and implement the *_impl() functions.

    */
    template<class ModelT = int> // dummy default to allow easy access to public typedefs
    class RansacModel
    {
  public:
      typedef std::pair<double, double> DPair;
      typedef std::vector<DPair> DVec;
      typedef DVec::const_iterator DVecIt;
      typedef std::vector<double> ModelParameters;


      /// fit a model and return its parameters
      ModelParameters rm_fit(const DVecIt& begin, const DVecIt& end) const
      {
        return static_cast<const ModelT*>(this)->rm_fit_impl(begin, end);
      }
      
      /**
        @brief Returns the R-squared of the data applied to the model (computed on-the-fly).
      
        Takes as input a standard vector of a standard pair of points in a 2D space.
        
        @param begin Iterator to first pair
        @param end Past-end iterator to last pair
        @return: R-squared value
      */
      double rm_rsq(const DVecIt& begin, const DVecIt& end) const
      {
        return static_cast<const ModelT*>(this)->rm_rsq_impl(begin, end);
      }

      /// calculates the residual sum of squares of the input points according to the model
      double rm_rss(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients) const
      {
        return static_cast<const ModelT*>(this)->rm_rss_impl(begin, end, coefficients);
      }

      /// calculates the squared residual of each input point vs. the model.
      /// All points with an error <= max_threshold, are accepted as inliers and returned.
      DVec rm_inliers(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold) const
      {
        return static_cast<const ModelT*>(this)->rm_inliers_impl(begin, end, coefficients, max_threshold);
      }
    };

  }


}
