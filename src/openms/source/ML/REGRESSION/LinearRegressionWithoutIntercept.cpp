// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

//#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/ML/REGRESSION/LinearRegressionWithoutIntercept.h>

namespace OpenMS::Math
{

    LinearRegressionWithoutIntercept::LinearRegressionWithoutIntercept() :
    sum_xx_(0), sum_xy_(0), n_(0)
    {
    }
    
    void LinearRegressionWithoutIntercept::addData(double x, double y)
    {
      sum_xx_ += x * x;
      sum_xy_ += x * y;
      
      ++n_;
    }
    
    void LinearRegressionWithoutIntercept::addData(std::vector<double>& x, std::vector<double>& y)
    {
      for (unsigned i = 0; i < x.size(); ++i)
      {
        addData(x[i], y[i]);
      }
    }

    /**
     * @brief returns the slope of the estimated regression line.
     */
    double LinearRegressionWithoutIntercept::getSlope() const
    {
      if (n_ < 2)
      {
        return std::numeric_limits<double>::quiet_NaN(); // not enough data
      }
      
      return sum_xy_ / sum_xx_;
    }

} //OpenMS //Math

