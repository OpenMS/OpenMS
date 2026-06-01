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
    @brief Natural cubic-spline interpolation of a 2D data set.

    Fits a piecewise cubic polynomial through the supplied (x, y)
    knots so that the resulting spline is twice continuously
    differentiable. The construction follows R.L. Burden, J.D. Faires,
    @e Numerical @e Analysis, 4th ed., PWS-Kent 1989,
    ISBN 0-53491-585-X, pp. 126-131.

    Construction is the expensive step; subsequent evaluations are
    roughly an order of magnitude (~50x) faster than construction in
    typical cases. After construction, the spline can be sampled at
    any @c x in the closed interval @c [x_first, x_last] via
    @ref eval, and its first three derivatives are available via
    @ref derivative / @ref derivatives.

    @ingroup Math
  */
  class OPENMS_DLLAPI CubicSpline2d
  {

    std::vector<double> a_; ///< Per-segment constant coefficient (degree 0); same value as the @c y at the segment's left knot.
    std::vector<double> b_; ///< Per-segment linear coefficient (degree 1).
    std::vector<double> c_; ///< Per-segment quadratic coefficient (degree 2); contains one extra trailing zero from the natural boundary condition.
    std::vector<double> d_; ///< Per-segment cubic coefficient (degree 3).
    std::vector<double> x_; ///< Sorted knot abscissae (length = number of segments + 1).

public:

    /**
      @brief Build the spline from parallel @p x / @p y vectors.

      @param[in] x Knot abscissae; must contain at least two entries
                   and be sorted in non-decreasing order.
      @param[in] y Knot ordinates; must have the same length as @p x.

      @throws Exception::IllegalArgument when @p x and @p y differ in
                                         length, when fewer than two
                                         knots are supplied, or when
                                         @p x is not non-decreasing.
    */
    CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y);

    /**
      @brief Build the spline from a map of (x, y) knots.

      The map is automatically sorted by @c x; passing an empty or
      single-entry map is an error.

      @param[in] m Knots as @c std::map<x, y>. Must contain at least two entries.

      @throws Exception::IllegalArgument when @p m contains fewer than
                                         two entries.
    */
    CubicSpline2d(const std::map<double, double>& m);

    /**
      @brief Evaluate the spline at @p x.

      @param[in] x Position; must lie in the closed interval @c [first_knot, last_knot].
      @return Interpolated value at @p x.

      @throws Exception::IllegalArgument when @p x lies outside the
                                         knot range.
    */
    double eval(double x) const;

    /**
      @brief First derivative of the spline at @p x.

      Equivalent to @c derivatives(x, @c 1).

      @param[in] x Position; must lie in the closed interval @c [first_knot, last_knot].
      @return First derivative at @p x.

      @throws Exception::IllegalArgument when @p x lies outside the
                                         knot range.
    */
    double derivative(double x) const;

    /**
      @brief Derivative of the spline at @p x.

      @param[in] x     Position; must lie in the closed interval
                       @c [first_knot, last_knot].
      @param[in] order Derivative order. Cubic splines provide
                       meaningful first, second, and third derivatives;
                       this method accepts @c order in @c {1, 2, 3}.
                       Higher orders are zero everywhere and are
                       rejected with an exception.
      @return @p order -th derivative of the spline at @p x.

      @throws Exception::IllegalArgument when @p x lies outside the
                                         knot range, or when @p order
                                         is not @c 1, @c 2, or @c 3.
    */
    double derivatives(double x, unsigned order) const;

private:

    /**
      @brief Initialise the spline coefficients from the supplied knots.

      Called from both public constructors after their argument
      validation; assumes the caller has already verified that the
      inputs are well-formed.

      @param[in] x Knot abscissae.
      @param[in] y Knot ordinates.
    */
    void init_(const std::vector<double>& x, const std::vector<double>& y);

  };

}

