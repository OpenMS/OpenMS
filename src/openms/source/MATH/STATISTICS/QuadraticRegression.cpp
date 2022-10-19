// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/MATH/STATISTICS/QuadraticRegression.h>

#include <Mathematics/Vector2.h>
#include <Mathematics/Vector3.h>
#include <Mathematics/Matrix3x3.h>
#include <Mathematics/LinearSystem.h>

namespace OpenMS::Math
{
  // Note:x, y must be of same size
  double computeChiSquareWeighted(
    std::vector<double>::const_iterator x_begin, 
    const std::vector<double>::const_iterator& x_end, 
    std::vector<double>::const_iterator y_begin,
    std::vector<double>::const_iterator w_begin,
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


  QuadraticRegression::QuadraticRegression() :
    a_(0), b_(0), c_(0), chi_squared_(0) {}

  double QuadraticRegression::eval(double x) const {return a_ + b_*x + c_*x*x;}

  double QuadraticRegression::eval(double A, double B, double C, double x) {return A + B*x + C*x*x;}

  double QuadraticRegression::getA() const {return a_;}
  double QuadraticRegression::getB() const {return b_;}
  double QuadraticRegression::getC() const {return c_;}
  double QuadraticRegression::getChiSquared() const {return chi_squared_;}


  void QuadraticRegression::computeRegression(std::vector<double>::const_iterator x_begin, 
    std::vector<double>::const_iterator x_end, 
    std::vector<double>::const_iterator y_begin)
  {
    std::vector<double> weights(std::distance(x_begin, x_end), 1);
    computeRegressionWeighted(x_begin, x_end, y_begin, weights.begin());
  }

  void QuadraticRegression::computeRegressionWeighted(
    std::vector<double>::const_iterator x_begin, 
    std::vector<double>::const_iterator x_end, 
    std::vector<double>::const_iterator y_begin, 
    std::vector<double>::const_iterator w_begin)
  {
    // Compute the linear fit of a quadratic function.
    // Get the coefficients for y = w_1*a +w_2*b*x + w_3*c*x^2.

    std::vector<gte::Vector2<double>> points;
    for(std::vector<double>::const_iterator xIter = x_begin, yIter = y_begin; xIter!=x_end; ++xIter, ++yIter)
    {
      points.emplace_back(std::initializer_list<double>{*xIter, *yIter});
    }

    // Compute sums for linear system. copy&paste from GeometricTools
    // and modified to allow quadratic functions
    int numPoints = static_cast<Int>(points.size());
    double sumX = 0, sumXX = 0, sumXXX = 0, sumXXXX = 0;
    double sumY = 0, sumXY = 0, sumXXY = 0;
    double sumW = 0;

    auto wIter = w_begin;
    for (int i = 0; i < numPoints; ++i, ++wIter)
    {
      double x = points[i][0];
      double y = points[i][1];
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
    gte::Matrix3x3<double> A
    {
      sumW, sumX, sumXX,
      sumX, sumXX, sumXXX,
      sumXX, sumXXX, sumXXXX
    };

    gte::Vector3<double> B
    {
      sumY,
      sumXY,
      sumXXY
    };

    gte::Vector3<double> X;

    bool nonsingular = gte::LinearSystem<double>().Solve(A, B, X);
    if (nonsingular)
    {
      a_ = X[0];
      b_ = X[1];
      c_ = X[2];
      chi_squared_ = computeChiSquareWeighted(x_begin, x_end, y_begin, w_begin, a_, b_, c_);
    }
    else
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-QuadraticRegression", "Could not fit a linear model to the data");
    }
  }
} //OpenMS //Math
