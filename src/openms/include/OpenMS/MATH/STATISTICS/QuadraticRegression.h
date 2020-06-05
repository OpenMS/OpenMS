// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Christian Ehrlich, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/STATISTICS/RegressionUtils.h>

#include "Wm5Vector2.h"
#include "Wm5LinearSystem.h"
#include <iterator>

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
      template <typename Iterator>
      void computeRegression(Iterator x_begin, Iterator x_end, Iterator y_begin);

      /** compute the weighted quadratic regression over 2D points */
      template <typename Iterator>
      void computeRegressionWeighted(Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin);

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

    namespace
    {
      //x, y must be of same size
      template <typename Iterator>
      double computeChiSquareWeighted_(
        Iterator x_begin, const Iterator& x_end, Iterator y_begin, Iterator w_begin,
        const double a, const double b, const double c)
      {
        double chi_squared(0.0);
        for (; x_begin != x_end; ++x_begin, ++y_begin, ++w_begin)
        {
          const double& x = *x_begin;
          const double& y = *y_begin;
          const double& weight = *w_begin;
          chi_squared += weight * std::pow(y - a - b * x - c * x * x, 2);
        }

        return chi_squared;
      }

    } // anonymous namespace

    template <typename Iterator>
    void QuadraticRegression::computeRegression(Iterator x_begin, Iterator x_end, Iterator y_begin)
    {
      std::vector<double> weights(std::distance(x_begin, x_end), 1);
      computeRegressionWeighted<Iterator>(x_begin, x_end, y_begin, weights.begin());
    }

    template <typename Iterator>
    void QuadraticRegression::computeRegressionWeighted(
      Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)
    {
      // Compute the linear fit of a quadratic function.
      // Get the coefficients for y = w_1*a +w_2*b*x + w_3*c*x^2.
      std::vector<Wm5::Vector2d> points = iteratorRange2Wm5Vectors(x_begin, x_end, y_begin);
      // Compute sums for linear system. copy&paste from GeometricTools Wm5ApprLineFit2.cpp
      // and modified to allow quadratic functions
      int numPoints = static_cast<Int>(points.size());
      double sumX = 0, sumXX = 0, sumXXX = 0, sumXXXX = 0;
      double sumY = 0, sumXY = 0, sumXXY = 0;
      double sumW = 0;

      Iterator wIter = w_begin;
      for (int i = 0; i < numPoints; ++i, ++wIter)
      {

        double x = points[i].X();
        double y = points[i].Y();
        double weight = *wIter;

        sumX += weight * x;
        sumXX += weight * x * x;
        sumXXX += weight * x * x * x;
        sumXXXX += weight * x * x * x * x;

        sumY += weight * y;
        sumXY += weight * x * y;
        sumXXY += weight * x * x * y;

        sumW += weight;
      }
      //create matrices to solve Ax = B
      double A[3][3] =
      {
        {sumW, sumX, sumXX},
        {sumX, sumXX, sumXXX},
        {sumXX, sumXXX, sumXXXX}
      };
      double B[3] =
      {
        sumY,
        sumXY,
        sumXXY
      };
      double X[3] = {0, 0, 0};

      bool nonsingular = Wm5::LinearSystem<double>().Solve3(A, B, X);
      if (nonsingular)
      {
        a_ = X[0];
        b_ = X[1];
        c_ = X[2];
        chi_squared_ = computeChiSquareWeighted_(x_begin, x_end, y_begin, w_begin, a_, b_, c_);
      }
      else
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-QuadraticRegression", "Could not fit a linear model to the data");
      }
    }

  } //namespace
} //namespace

