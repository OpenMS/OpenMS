// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>

#include <cmath>
#include <vector>

namespace OpenMS
{
  namespace Math
  {
    /**
      @brief This class offers functions to perform least-squares fits to a straight line model, \f$ Y(c,x) = c_0 + c_1 x \f$.

      @ingroup Math
    */
    class OPENMS_DLLAPI LinearRegressionWithoutIntercept
    {
      public:

      /// Constructor
      LinearRegressionWithoutIntercept();
      
      /**
       * @brief adds an observation (x,y) to the regression data set.
       *
       * @param x    independent variable value
       * @param y    dependent variable value
       */
      void addData(double x, double y);
      
      /**
       * @brief adds observations (x,y) to the regression data set.
       *
       * @param x    vector of independent variable values
       * @param y    vector of dependent variable values
       */
      void addData(std::vector<double>& x, std::vector<double>& y);
      
      /**
       * @brief returns the slope of the estimated regression line.
       */
      double getSlope() const;

      private:
      /**
       * @brief total variation in x
       */
      double sum_xx_;
      
      /**
       * @brief sum of products
       */
      double sum_xy_;
      
      /**
       * @brief number of observations
       */
      int n_;
      
    };

  } // namespace Math
} // namespace OpenMS


