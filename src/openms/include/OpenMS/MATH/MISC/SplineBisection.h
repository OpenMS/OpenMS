// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <limits> 
#include <cmath> 

namespace OpenMS
{
  /**
   * @brief Uses bisection to find the maximum point of a spline
   *
   * Should work with BSpline2d and CubicSpline2d
   *
   */
  namespace Math
  {

    template <class T>
    void spline_bisection(const T & peak_spline, 
        double const left_neighbor_mz,
        double const right_neighbor_mz,
        double & max_peak_mz,
        double & max_peak_int,
        double const threshold = 1e-6)
    {
      // calculate maximum by evaluating the spline's 1st derivative
      // (bisection method)
      double lefthand = left_neighbor_mz;
      double righthand = right_neighbor_mz;

      bool lefthand_sign = true;
      double eps = std::numeric_limits<double>::epsilon();

      // bisection
      do
      {
        double mid = (lefthand + righthand) / 2.0;
        double midpoint_deriv_val = peak_spline.derivative(mid);

        // if deriv nearly zero then maximum already found
        if (!(std::fabs(midpoint_deriv_val) > eps))
        {
          break;
        }

        bool midpoint_sign = (midpoint_deriv_val < 0.0) ? false : true;

        if (lefthand_sign ^ midpoint_sign)
        {
          righthand = mid;
        }
        else
        {
          lefthand = mid;
        }
      }
      while (righthand - lefthand > threshold);

      max_peak_mz = (lefthand + righthand) / 2;
      max_peak_int = peak_spline.eval(max_peak_mz);
    }

  }
}
