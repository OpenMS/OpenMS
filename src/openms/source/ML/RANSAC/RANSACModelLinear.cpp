// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ML/RANSAC/RANSACModelLinear.h>

#include <OpenMS/ML/REGRESSION/LinearRegression.h>
#include <numeric>


namespace OpenMS::Math
{

    RansacModelLinear::ModelParameters RansacModelLinear::rm_fit_impl(const DVecIt& begin, const DVecIt& end)
    {
      std::vector<double> x, y;

      for (DVecIt it = begin; it != end; ++it)
      {
        x.push_back(it->first);
        y.push_back(it->second);
      }
      LinearRegression lin_reg;
      lin_reg.computeRegression(0.95, x.begin(), x.end(), y.begin(), false); // no goodness of fit computation
      ModelParameters p;
      p.push_back(lin_reg.getIntercept());
      p.push_back(lin_reg.getSlope());
      return p;
    }

    double RansacModelLinear::rm_rsq_impl(const DVecIt& begin, const DVecIt& end)
    {
      std::vector<double> x, y;

      for (DVecIt it = begin; it != end; ++it)
      {
        x.push_back(it->first);
        y.push_back(it->second);
      }

      LinearRegression lin_reg;
      lin_reg.computeRegression(0.95, x.begin(), x.end(), y.begin(), false);

      return lin_reg.getRSquared();
    }

    double RansacModelLinear::rm_rss_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients)
    {
      double rss = 0;

      for (DVecIt it = begin; it != end; ++it)
      {
        rss += pow(it->second - (coefficients[0] + ( coefficients[1] * it->first)), 2);
      }

      return rss;
    }

    RansacModelLinear::DVec RansacModelLinear::rm_inliers_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold)
    {
      DVec alsoinliers;
      //std::cerr << "\n\nRANSAC dists: ";
      for (DVecIt it = begin; it != end; ++it)
      {
        double dist = pow(it->second - (coefficients[0] + ( coefficients[1] * it->first)), 2);
        //std::cerr << dist << ", ";
        if (dist < max_threshold)
        {
          alsoinliers.push_back(*it);
        }
      }
      //std::cerr << "\n\n";

      return alsoinliers;
    }


} // OpenMS //Math
