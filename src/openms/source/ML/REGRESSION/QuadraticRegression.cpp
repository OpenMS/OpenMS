// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Christian Ehrlich, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ML/REGRESSION/QuadraticRegression.h>

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
