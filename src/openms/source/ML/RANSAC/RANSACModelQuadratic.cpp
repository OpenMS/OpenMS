// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ML/RANSAC/RANSACModelQuadratic.h>

#include <OpenMS/ML/REGRESSION/QuadraticRegression.h>
#include <numeric>


namespace OpenMS::Math
{

    // Quadratic regression for RANSAC
    RansacModelQuadratic::ModelParameters RansacModelQuadratic::rm_fit_impl(const DVecIt& begin, const DVecIt& end)
    {
      std::vector<double> x, y;

      for (DVecIt it = begin; it != end; ++it)
      {
        x.push_back(it->first);
        y.push_back(it->second);
      }
      QuadraticRegression quad_reg;
      quad_reg.computeRegression(x.begin(), x.end(), y.begin());
      ModelParameters p;
      p.push_back(quad_reg.getA());
      p.push_back(quad_reg.getB());
      p.push_back(quad_reg.getC());
      return p;
    }

    double RansacModelQuadratic::rm_rsq_impl(const DVecIt& begin, const DVecIt& end)
    {
      std::vector<double> x, y;

      for (DVecIt it = begin; it != end; ++it)
      {
        x.push_back(it->first);
        y.push_back(it->second);
      }

      QuadraticRegression quad_reg;
      quad_reg.computeRegression(x.begin(), x.end(), y.begin());

      return quad_reg.getChiSquared();
    }

    double RansacModelQuadratic::rm_rss_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients)
    {
      double rss = 0;

      for (DVecIt it = begin; it != end; ++it)
      {
        const double value_model = QuadraticRegression::eval(coefficients[0], coefficients[1], coefficients[2], it->first);
        const double diff = it->second - value_model;
        rss += diff*diff;
      }

      return rss;
    }

    RansacModelQuadratic::DVec RansacModelQuadratic::rm_inliers_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold)
    {
      DVec alsoinliers;

      for (DVecIt it = begin; it != end; ++it)
      {
        const double value_model = QuadraticRegression::eval(coefficients[0], coefficients[1], coefficients[2], it->first);
        const double diff = it->second - value_model;
        if (diff * diff < max_threshold)
        {
          alsoinliers.push_back(*it);
        }
      }

      return alsoinliers;
    }


} // OpenMS //Math
