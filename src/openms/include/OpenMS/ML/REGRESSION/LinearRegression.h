// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

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
      virtual ~LinearRegression() = default;

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

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      void computeRegression(double confidence_interval_P,
        std::vector<double>::const_iterator x_begin,
        std::vector<double>::const_iterator x_end,
        std::vector<double>::const_iterator y_begin,
        bool compute_goodness = true);

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

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      void computeRegressionWeighted(double confidence_interval_P, 
        std::vector<double>::const_iterator x_begin, 
        std::vector<double>::const_iterator x_end, 
        std::vector<double>::const_iterator y_begin, 
        std::vector<double>::const_iterator w_begin, 
        bool compute_goodness = true);

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

      /// given x compute y = slope * x + intercept
      static inline double computePointY(double x, double slope, double intercept)
      {
        return slope * x + intercept;
      }

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
      void computeGoodness_(const std::vector<double>& X, const std::vector<double>& Y, double confidence_interval_P);

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
  } // namespace Math
} // namespace OpenMS


