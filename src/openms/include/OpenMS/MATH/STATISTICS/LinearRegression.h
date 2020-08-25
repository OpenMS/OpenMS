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
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/MATH/STATISTICS/RegressionUtils.h>

#include "Wm5Vector2.h"
#include "Wm5ApprLineFit2.h"
#include "Wm5LinearSystem.h"

#include <cmath>
#include <vector>

namespace OpenMS
{
  namespace Math
  {
    /**
      @brief This class offers functions to perform least-squares fits to a straight line model, \f$ Y(c,x) = c_0 + c_1 x \f$.

      Next to the intercept with the y-axis and the slope of the fitted line, this class computes the:
        - squared Pearson coefficient
        - value of the t-distribution
        - standard deviation of the residuals
        - standard error of the slope
        - intercept with the x-axis (useful for additive series experiments)
        - lower border of confidence interval
        - higher border of confidence interval
        - chi squared value
        - x mean

      @ingroup Math
    */
    class OPENMS_DLLAPI LinearRegression
    {
public:

      /// Constructor
      LinearRegression() :
        intercept_(0),
        slope_(0),
        x_intercept_(0),
        lower_(0),
        upper_(0),
        t_star_(0),
        r_squared_(0),
        stand_dev_residuals_(0),
        mean_residuals_(0),
        stand_error_slope_(0),
        chi_squared_(0),
        rsd_(0)
      {
      }

      /// Destructor
      virtual ~LinearRegression()
      {
      }

      /**
          @brief This function computes the best-fit linear regression coefficients \f$ (c_0,c_1) \f$
          of the model \f$ Y = c_0 + c_1 X \f$ for the dataset \f$ (x, y) \f$.

          The values in x-dimension of the dataset \f$ (x,y) \f$ are given by the iterator range [x_begin,x_end)
          and the corresponding y-values start at position y_begin.


          For a  "x %" Confidence Interval use confidence_interval_P = x/100.
          For example the 95% Confidence Interval is supposed to be an interval that has a 95% chance of
          containing the true value of the parameter.
          
          @param confidence_interval_P Value between 0-1 to determine lower and upper CI borders.
          @param x_begin Begin iterator of x values
          @param x_end End iterator of x values
          @param y_begin Begin iterator of y values (same length as x)
          @param compute_goodness Compute meta stats about the fit. If this is not done, none of the members (except slope and intercept) are meaningful.
          @return If an error occurred during the fit.

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      template <typename Iterator>
      void computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, bool compute_goodness = true);

      /**
          @brief This function computes the best-fit linear regression coefficients \f$ (c_0,c_1) \f$
          of the model \f$ Y = c_0 + c_1 X \f$ for the weighted dataset \f$ (x, y) \f$.

          The values in x-dimension of the dataset \f$ (x, y) \f$ are given by the iterator range [x_begin,x_end)
          and the corresponding y-values start at position y_begin. They will be weighted by the
          values starting at w_begin.

          For a  "x %" Confidence Interval use confidence_interval_P = x/100.
          For example the 95% Confidence Interval is supposed to be an interval that has a 95% chance of
          containing the true value of the parameter.

          @param confidence_interval_P Value between 0-1 to determine lower and upper CI borders.
          @param x_begin Begin iterator of x values
          @param x_end End iterator of x values
          @param y_begin Begin iterator of y values (same length as x)
          @param w_begin Begin iterator of weight values (same length as x)
          @param compute_goodness Compute meta stats about the fit. If this is not done, none of the members (except slope and intercept) are meaningful.
          @return If an error occurred during the fit.

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      template <typename Iterator>
      void computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin, bool compute_goodness = true);

      /// Non-mutable access to the y-intercept of the straight line
      double getIntercept() const;
      /// Non-mutable access to the slope of the straight line
      double getSlope() const;
      /// Non-mutable access to the x-intercept of the straight line
      double getXIntercept() const;
      /// Non-mutable access to the lower border of confidence interval
      double getLower() const;
      /// Non-mutable access to the upper border of confidence interval
      double getUpper() const;
      /// Non-mutable access to the value of the t-distribution
      double getTValue() const;
      /// Non-mutable access to the squared Pearson coefficient
      double getRSquared() const;
      /// Non-mutable access to the standard deviation of the residuals
      double getStandDevRes() const;
      /// Non-mutable access to the residual mean
      double getMeanRes() const;
      /// Non-mutable access to the standard error of the slope
      double getStandErrSlope() const;
      /// Non-mutable access to the chi squared value
      double getChiSquared() const;
      /// Non-mutable access to relative standard deviation
      double getRSD() const;

protected:

      /// The intercept of the fitted line with the y-axis
      double intercept_;
      /// The slope of the fitted line
      double slope_;
      /// The intercept of the fitted line with the x-axis
      double x_intercept_;
      /// The lower bound of the confidence interval
      double lower_;
      /// The upper bound of the confidence interval
      double upper_;
      /// The value of the t-statistic
      double t_star_;
      /// The squared correlation coefficient (Pearson)
      double r_squared_;
      /// The standard deviation of the residuals
      double stand_dev_residuals_;
      /// Mean of residuals
      double mean_residuals_;
      /// The standard error of the slope
      double stand_error_slope_;
      /// The value of the Chi Squared statistic
      double chi_squared_;
      /// the relative standard deviation
      double rsd_;


      /// Computes the goodness of the fitted regression line
      void computeGoodness_(const std::vector<Wm5::Vector2d>& points, double confidence_interval_P);

      /// Compute the chi squared of a linear fit
      template <typename Iterator>
      double computeChiSquare(Iterator x_begin, Iterator x_end, Iterator y_begin, double slope, double intercept);

      /// Compute the chi squared of a weighted linear fit
      template <typename Iterator>
      double computeWeightedChiSquare(Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin, double slope, double intercept);

private:

      /// Not implemented
      LinearRegression(const LinearRegression& arg);
      /// Not implemented
      LinearRegression& operator=(const LinearRegression& arg);

    }; //class


    namespace
    {
      //given x compute y = slope * x + intercept
      double computePointY(double x, double slope, double intercept)
      {
        return slope * x + intercept;
      }

    } //namespace

    //x, y, w must be of same size
    template <typename Iterator>
    double LinearRegression::computeChiSquare(Iterator x_begin, Iterator x_end, Iterator y_begin, double slope, double intercept)
    {
      double chi_squared = 0.0;
      Iterator xIter = x_begin;
      Iterator yIter = y_begin;
      for (; xIter != x_end; ++xIter, ++yIter)
      {
        chi_squared += std::pow(*yIter - computePointY(*xIter, slope, intercept), 2);
      }

      return chi_squared;
    }

    //x, y, w must be of same size
    template <typename Iterator>
    double LinearRegression::computeWeightedChiSquare(Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin, double slope, double intercept)
    {
      double chi_squared = 0.0;
      Iterator xIter = x_begin;
      Iterator yIter = y_begin;
      Iterator wIter = w_begin;
      for (; xIter != x_end; ++xIter, ++yIter, ++wIter)
      {
        chi_squared += *wIter * std::pow(*yIter - computePointY(*xIter, slope, intercept), 2);
      }

      return chi_squared;
    }

    template <typename Iterator>
    void LinearRegression::computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, bool compute_goodness)
    {
      std::vector<Wm5::Vector2d> points = iteratorRange2Wm5Vectors(x_begin, x_end, y_begin);

      // Compute the unweighted linear fit.
      // Get the intercept and the slope of the regression Y_hat=intercept_+slope_*X
      // and the value of Chi squared (sum( (y - evel(x))^2)
      bool pass = Wm5::HeightLineFit2<double>(static_cast<int>(points.size()), &points.front(), slope_, intercept_);
      chi_squared_ = computeChiSquare(x_begin, x_end, y_begin, slope_, intercept_);

      if (pass)
      {
        if (compute_goodness && points.size() > 2) computeGoodness_(points, confidence_interval_P);
      }
      else
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "UnableToFit-LinearRegression", String("Could not fit a linear model to the data (") + points.size() + " points).");
      }
    }

    template <typename Iterator>
    void LinearRegression::computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin, bool compute_goodness)
    {
      // Compute the weighted linear fit.
      // Get the intercept and the slope of the regression Y_hat=intercept_+slope_*X
      // and the value of Chi squared, the covariances of the intercept and the slope
      std::vector<Wm5::Vector2d> points = iteratorRange2Wm5Vectors(x_begin, x_end, y_begin);
      // Compute sums for linear system. copy&paste from GeometricTools Wm5ApprLineFit2.cpp
      // and modified to allow weights
      int numPoints = static_cast<int>(points.size());
      double sumX = 0, sumY = 0;
      double sumXX = 0, sumXY = 0;
      double sumW = 0;
      Iterator wIter = w_begin;

      for (int i = 0; i < numPoints; ++i, ++wIter)
      {
        sumX += (*wIter) * points[i].X();
        sumY += (*wIter) * points[i].Y();
        sumXX += (*wIter) * points[i].X() * points[i].X();
        sumXY += (*wIter) * points[i].X() * points[i].Y();
        sumW += (*wIter);
      }
      //create matrices to solve Ax = B
      double A[2][2] =
      {
        {sumXX, sumX},
        {sumX, sumW}
      };
      double B[2] =
      {
        sumXY,
        sumY
      };
      double X[2];

      bool nonsingular = Wm5::LinearSystem<double>().Solve2(A, B, X);
      if (nonsingular)
      {
        slope_ = X[0];
        intercept_ = X[1];
      }
      chi_squared_ = computeWeightedChiSquare(x_begin, x_end, y_begin, w_begin, slope_, intercept_);

      if (nonsingular)
      {
        if (compute_goodness && points.size() > 2) computeGoodness_(points, confidence_interval_P);
      }
      else
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "UnableToFit-LinearRegression", "Could not fit a linear model to the data");
      }
    }

  } // namespace Math
} // namespace OpenMS


