// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Christian Ehrlich, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iterator>
#include <vector>

namespace OpenMS
{
  namespace Math
  {
    /**
      @brief Estimates model parameters for a quadratic equation.

      Fits the coefficients of
       @c y = a + b*x + c*x^2
      via the normal equations of a (weighted) least-squares problem.
      Both an unweighted (@ref computeRegression) and a weighted
      (@ref computeRegressionWeighted) overload are provided.

      Once a fit is computed the coefficients and the chi-squared sum
      can be queried via @ref getA / @ref getB / @ref getC / @ref getChiSquared,
      and the resulting polynomial can be evaluated via @ref eval. A
      static overload of @ref eval is available for evaluating an
      externally-supplied @c (A,B,C) triple.

      Before the first successful fit the coefficients and the
      chi-squared sum read back as @c 0.

      @ingroup Math
    */
    class OPENMS_DLLAPI QuadraticRegression
    {
public:
      /**
        @brief Default constructor.

        Leaves all coefficients and the chi-squared sum at @c 0; the
        instance produces meaningful results only after the first
        successful call to @ref computeRegression or
        @ref computeRegressionWeighted.
      */
      QuadraticRegression();

      /**
        @brief Fit @c y = a + b*x + c*x^2 to unweighted 2D points.

        Equivalent to calling @ref computeRegressionWeighted with a
        uniform weight of @c 1 for every point. On success the
        coefficients and @ref getChiSquared are updated.

        @param[in] x_begin Iterator to the first @c x value.
        @param[in] x_end   Past-the-end iterator for the @c x range.
        @param[in] y_begin Iterator to the first @c y value; the @c y range must contain at least @c std::distance(x_begin, x_end) elements.

        @throws Exception::UnableToFit if the resulting linear system is singular (e.g. fewer than three distinct @c x values).
      */
      void computeRegression(
        std::vector<double>::const_iterator x_begin,
        std::vector<double>::const_iterator x_end,
        std::vector<double>::const_iterator y_begin);

      /**
        @brief Fit @c y = a + b*x + c*x^2 to weighted 2D points.

        Solves the weighted normal equations for the quadratic model. On
        success the coefficients and @ref getChiSquared are updated.

        @param[in] x_begin Iterator to the first @c x value.
        @param[in] x_end   Past-the-end iterator for the @c x range.
        @param[in] y_begin Iterator to the first @c y value; the @c y range must contain at least @c std::distance(x_begin, x_end) elements.
        @param[in] w_begin Iterator to the first weight; the weight range must contain at least @c std::distance(x_begin, x_end) elements.

        @throws Exception::UnableToFit if the resulting linear system is singular (e.g. fewer than three distinct @c x values).
      */
      void computeRegressionWeighted(
        std::vector<double>::const_iterator x_begin,
        std::vector<double>::const_iterator x_end,
        std::vector<double>::const_iterator y_begin,
        std::vector<double>::const_iterator w_begin);

      /**
        @brief Evaluate the fitted polynomial at @p x.

        @param[in] x Input value.
        @return @c a + b*x + c*x^2 using the currently stored coefficients.
      */
      double eval(double x) const;

      /**
        @brief Evaluate a quadratic polynomial at @p x using externally-supplied coefficients.

        @param[in] A Constant term.
        @param[in] B Linear coefficient.
        @param[in] C Quadratic coefficient.
        @param[in] x Input value.
        @return @c A + B*x + C*x*x.
      */
      static double eval(double A, double B, double C, double x);

      /**
        @brief Constant term of the fitted polynomial.

        @return @c a from the last successful fit (or @c 0 if no fit has been computed yet).
      */
      double getA() const;

      /**
        @brief Linear coefficient of the fitted polynomial.

        @return @c b from the last successful fit (or @c 0 if no fit has been computed yet).
      */
      double getB() const;

      /**
        @brief Quadratic coefficient of the fitted polynomial.

        @return @c c from the last successful fit (or @c 0 if no fit has been computed yet).
      */
      double getC() const;

      /**
        @brief Weighted residual sum of squares of the last successful fit.

        @return @c sum(weight * (y - a - b*x - c*x*x)^2) over the last
                fit's input points (or @c 0 if no fit has been computed
                yet).
      */
      double getChiSquared() const;

protected:
      double a_;           ///< Constant term of the last successful fit.
      double b_;           ///< Linear coefficient of the last successful fit.
      double c_;           ///< Quadratic coefficient of the last successful fit.
      double chi_squared_; ///< Weighted residual sum of the last successful fit.
    }; //class

  } //namespace
} //namespace

