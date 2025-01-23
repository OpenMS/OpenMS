// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>

#include <vector>
#include <map>

namespace OpenMS
{
  /**
    @brief cubic spline interpolation
    as described in R.L. Burden, J.D. Faires, Numerical Analysis, 4th ed.
    PWS-Kent, 1989, ISBN 0-53491-585-X, pp. 126-131.
    
    Construction of the spline takes by far the most time. Evaluating it is rather fast 
    (one evaluation is about 50x faster than construction -- depending on number of points etc.).
   
   */
  class OPENMS_DLLAPI CubicSpline2d
  {

    std::vector<double> a_; ///< constant spline coefficients
    std::vector<double> b_; ///< linear spline coefficients
    std::vector<double> c_; ///< quadratic spline coefficients
    std::vector<double> d_; ///< cubic spline coefficients
    std::vector<double> x_; ///< knots

public:

    /**
     * @brief constructor of spline interpolation
     *
     * The coordinates must match by index. Both vectors must be
     * the same size and sorted in x. Sortedness in x is required
     * for @see SplinePackage.
     *
     * @param x x-coordinates of input data points (knots)
     * @param y y-coordinates of input data points
     */
    CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y);

    /**
     * @brief constructor of spline interpolation
     *
     * @param m (x,y) coordinates of input data points
     */
    CubicSpline2d(const std::map<double, double>& m);

    /**
     * @brief evaluates the spline at position x
     *
     * @param x x-position
     */
    double eval(double x) const;

    /**
     * @brief evaluates first derivative of spline at position x
     *
     * @param x x-position
     */
    double derivative(double x) const;

    /**
     * @brief evaluates derivative of spline at position x
     *
     * @param x x-position
     * @param order order of the derivative
     * Only order 1 or 2 make sense for cubic splines.
     */
    double derivatives(double x, unsigned order) const;

private:

    /**
     * @brief initialize the spline
     *
     * @param x x-coordinates of input data points (knots)
     * @param y y-coordinates of input data points
     */
    void init_(const std::vector<double>& x, const std::vector<double>& y);

  };

}

