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
    /*
      @brief Estimates model parameters for a quadratic equation

      The quadratic equation is of the form 
       y = a + b*x + c*x^2

      Weighted inputs are supported. 

    */
    class OPENMS_DLLAPI QuadraticRegression
    {
public:
      QuadraticRegression();

      /** compute the quadratic regression over 2D points */
      void computeRegression(
        std::vector<double>::const_iterator x_begin, 
        std::vector<double>::const_iterator x_end, 
        std::vector<double>::const_iterator y_begin);

      /** compute the weighted quadratic regression over 2D points */
      void computeRegressionWeighted(
        std::vector<double>::const_iterator x_begin, 
        std::vector<double>::const_iterator x_end,
        std::vector<double>::const_iterator y_begin, 
        std::vector<double>::const_iterator w_begin);

      /** evaluate the quadratic function */
      double eval(double x) const;

      /** evaluate using external coefficients */
      static double eval(double A, double B, double C, double x);

      double getA() const; /// A = the intercept
      double getB() const; /// B*X
      double getC() const; /// C*X^2
      double getChiSquared() const;

protected:
      double a_;
      double b_;
      double c_;
      double chi_squared_;
    }; //class

  } //namespace
} //namespace

