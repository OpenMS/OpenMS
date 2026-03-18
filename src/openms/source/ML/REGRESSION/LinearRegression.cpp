// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/ML/REGRESSION/LinearRegression.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <Mathematics/Vector2.h>
#include <Mathematics/ApprHeightLine2.h>
#include <Mathematics/LinearSystem.h>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/distributions.hpp>

#include <iostream>

using boost::math::detail::inverse_students_t;

namespace OpenMS::Math
{
    static void vector2ToStdVec_(const std::vector<gte::Vector2<double>>& points, std::vector<double>& Xout, std::vector<double>& Yout){
      unsigned N = static_cast<unsigned>(points.size());
      Xout.clear(); Xout.reserve(N);
      Yout.clear(); Yout.reserve(N);
      for (unsigned i = 0; i < N; ++i)
      {
        Xout.push_back(points[i][0]);
        Yout.push_back(points[i][1]);
      }
    }

    double LinearRegression::getIntercept() const
    {
      return intercept_;
    }

    double LinearRegression::getSlope() const
    {
      return slope_;
    }

    double LinearRegression::getXIntercept() const
    {
      return x_intercept_;
    }

    double LinearRegression::getLower() const
    {
      return lower_;
    }

    double LinearRegression::getUpper() const
    {
      return upper_;
    }

    double LinearRegression::getTValue() const
    {
      return t_star_;
    }

    double LinearRegression::getRSquared() const
    {
      return r_squared_;
    }

    double LinearRegression::getStandDevRes() const
    {
      return stand_dev_residuals_;
    }

    double LinearRegression::getMeanRes() const
    {
      return mean_residuals_;
    }

    double LinearRegression::getStandErrSlope() const
    {
      return stand_error_slope_;
    }

    double LinearRegression::getChiSquared() const
    {
      return chi_squared_;
    }

    double LinearRegression::getRSD() const
    {
      return rsd_;
    }

    void LinearRegression::computeGoodness_(const std::vector<double>& X, const std::vector<double>& Y, double confidence_interval_P)
    {
      OPENMS_PRECONDITION(static_cast<unsigned>(X.size() == Y.size()), 
          "Fitted X and Y have different lengths.");
      OPENMS_PRECONDITION(static_cast<unsigned>(X.size()) > 2, 
          "Cannot compute goodness of fit for regression with less than 3 data points");
      // specifically, boost throws an exception for a t-distribution with zero df

      Size N = X.size();

      // Mean of abscissa and ordinate values
      double x_mean = Math::mean(X.begin(), X.end());
      double y_mean = Math::mean(Y.begin(), Y.end());

      // Variance and Covariances
      double var_X = Math::variance(X.begin(), X.end(), x_mean);
      double var_Y = Math::variance(Y.begin(), Y.end(), y_mean);
      double cov_XY = Math::covariance(X.begin(), X.end(), Y.begin(), Y.end());

      // S_xx
      double s_XX = var_X * (N-1);
      /*for (unsigned i = 0; i < N; ++i)
      {
        double d = (X[i] - x_mean);
        s_XX += d * d;
      }*/

      // Compute the squared Pearson coefficient
      r_squared_ = (cov_XY * cov_XY) / (var_X * var_Y);

      // The standard deviation of the residuals
      double sum = 0;
      for (unsigned i = 0; i < N; ++i)
      {
        double x_i = fabs(Y[i] - (intercept_ + slope_ * X[i]));
        sum += x_i;
      }
      mean_residuals_       = sum / N;
      stand_dev_residuals_ = sqrt((chi_squared_ - (sum * sum) / N) / (N - 1));

      // The Standard error of the slope
      stand_error_slope_ = stand_dev_residuals_ / sqrt(s_XX);

      // and the intersection of Y_hat with the x-axis
      x_intercept_ = -(intercept_ / slope_);

      double P = 1 - (1 - confidence_interval_P) / 2;
      boost::math::students_t tdist(N - 2);
      t_star_ = boost::math::quantile(tdist, P);

      //Compute the asymmetric 95% confidence interval of around the X-intercept
      double g = (t_star_ / (slope_ / stand_error_slope_));
      g *= g;
      double left = (x_intercept_ - x_mean) * g;
      double bottom = 1 - g;
      double d = (x_intercept_ - x_mean);
      double right = t_star_ * (stand_dev_residuals_ / slope_) * sqrt((d * d) / s_XX + (bottom / N));

      // Confidence interval lower_ <= X_intercept <= upper_
      lower_ = x_intercept_ + (left + right) / bottom;
      upper_ = x_intercept_ + (left - right) / bottom;

      if (lower_ > upper_)
      {
        std::swap(lower_, upper_);
      }

      double tmp = 0;
      for (unsigned i = 0; i < N; ++i)
      {
        tmp += (X[i] - x_mean) * (X[i] - x_mean);
      }

      //            cout << "100.0 / abs( x_intercept_ ) " << (100.0 / fabs( x_intercept_ )) << endl;
      //            cout << "tmp : " << tmp << endl;
      //            cout << "slope_ " << slope_ << endl;
      //            cout << "y_mean " << y_mean << endl;
      //            cout << "N " << N << endl;
      //            cout << "stand_dev_residuals_ " << stand_dev_residuals_ << endl;
      //            cout << " (1.0/ (double) N)  " <<  (1.0/ (double) N)  << endl;
      //            cout << "sx hat " << (stand_dev_residuals_ / slope_) * sqrt(  (1.0/ (double) N) * (y_mean / (slope_ * slope_ * tmp ) ) ) << endl;

      // compute relative standard deviation (non-standard formula, taken from Mayr et al. (2006) )
      rsd_ = (100.0 / fabs(x_intercept_)) * (stand_dev_residuals_ / slope_) * sqrt((1.0 / (double) N) * (y_mean / (slope_ * slope_ * tmp)));

      if (rsd_ < 0.0)
      {
        std::cout << "rsd < 0.0 " << std::endl;
        std::cout << "Intercept                                  " << intercept_
                  << "\nSlope                                    " << slope_
                  << "\nSquared pearson coefficient              " << r_squared_
                  << "\nValue of the t-distribution              " << t_star_
                  << "\nStandard deviation of the residuals      " << stand_dev_residuals_
                  << "\nStandard error of the slope              " << stand_error_slope_
                  << "\nThe X intercept                          " << x_intercept_
                  << "\nThe lower border of confidence interval  " << lower_
                  << "\nThe higher border of confidence interval " << upper_
                  << "\nChi squared value                        " << chi_squared_
                  << "\nx mean                                   " << x_mean
                  << "\nstand_error_slope/slope_                 " << (stand_dev_residuals_ / slope_)
                  << "\nCoefficient of Variation                 " << (stand_dev_residuals_ / slope_) / x_mean * 100  << std::endl
                  << "========================================="
                  << std::endl;
      }
    }

    void LinearRegression::computeRegression(double confidence_interval_P,
        std::vector<double>::const_iterator x_begin,
        std::vector<double>::const_iterator x_end,
        std::vector<double>::const_iterator y_begin,
        bool compute_goodness)
    {
      std::vector<gte::Vector2<double>> points;
      for(std::vector<double>::const_iterator xIter = x_begin, yIter = y_begin; xIter!=x_end; ++xIter, ++yIter)
      {
        points.emplace_back(std::initializer_list<double>{*xIter, *yIter});
      }

      // Compute the unweighted linear fit.
      // Get the intercept and the slope of the regression Y_hat=intercept_+slope_*X
      // and the value of Chi squared (sum( (y - evel(x))^2)
      auto line = gte::ApprHeightLine2<double>();
      bool pass = line.Fit(static_cast<int>(points.size()), &points.front());
      slope_ = line.GetParameters().second[0];
      intercept_ = -slope_ * line.GetParameters().first[0] + line.GetParameters().first[1];
      chi_squared_ = computeChiSquare(x_begin, x_end, y_begin, slope_, intercept_);

      if (!pass)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-LinearRegression", String("Could not fit a linear model to the data (") + points.size() + " points).");
      }

      if (compute_goodness && points.size() > 2)
      {
        std::vector<double> X,Y;
        vector2ToStdVec_(points, X, Y);
        computeGoodness_(X, Y , confidence_interval_P);
      }
    }

    void LinearRegression::computeRegressionWeighted(double confidence_interval_P, 
        std::vector<double>::const_iterator x_begin, 
        std::vector<double>::const_iterator x_end, 
        std::vector<double>::const_iterator y_begin, 
        std::vector<double>::const_iterator w_begin, 
        bool compute_goodness)    
     {
      // Compute the weighted linear fit.
      // Get the intercept and the slope of the regression Y_hat=intercept_+slope_*X
      // and the value of Chi squared, the covariances of the intercept and the slope
      std::vector<gte::Vector2<double>> points;
      for(std::vector<double>::const_iterator xIter = x_begin, yIter = y_begin; xIter!=x_end; ++xIter, ++yIter)
      {
        points.emplace_back(std::initializer_list<double>{*xIter, *yIter});
      }

      // Compute sums for linear system. copy&paste from GeometricToolsEngine (gte) ApprHeightLine2.h
      // and modified to allow weights
      int numPoints = static_cast<int>(points.size());
      double sumX = 0, sumY = 0;
      double sumXX = 0, sumXY = 0;
      double sumW = 0;
      auto wIter = w_begin;

      for (int i = 0; i < numPoints; ++i, ++wIter)
      {
        sumX += (*wIter) * points[i][0];
        sumY += (*wIter) * points[i][1];
        sumXX += (*wIter) * points[i][0] * points[i][0];
        sumXY += (*wIter) * points[i][0] * points[i][1];
        sumW += (*wIter);
      }
      //create matrices to solve Ax = B
      gte::Matrix2x2<double> A
      {
        sumXX, sumX,
        sumX, sumW
      };

      gte::Vector2<double> B
      {
        sumXY,
        sumY
      };

      gte::Vector2<double> X;

      bool nonsingular = gte::LinearSystem<double>().Solve(A, B, X);
      if (nonsingular)
      {
        slope_ = X[0];
        intercept_ = X[1];
      }
      chi_squared_ = computeWeightedChiSquare(x_begin, x_end, y_begin, w_begin, slope_, intercept_);

      if (nonsingular)
      {
        if (compute_goodness && points.size() > 2)
        {
          std::vector<double> X,Y;
          vector2ToStdVec_(points, X, Y);
          computeGoodness_(X, Y, confidence_interval_P);
        }
      }
      else
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "UnableToFit-LinearRegression", "Could not fit a linear model to the data");
      }
    }
} // OpenMS //Math

