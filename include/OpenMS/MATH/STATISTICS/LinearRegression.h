// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_LINEARREGRESSION_H
#define OPENMS_MATH_STATISTICS_LINEARREGRESSION_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>
#include <vector>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>

namespace OpenMS
{
  namespace Math
  {
    /**
			@brief This class offers functions to perform least-squares fits to a straight line model, \f$ Y(c,x) = c_0 + c_1 x \f$.

			It capsulates the GSL methods for a weighted and an unweighted linear regression.

			Next to the intercept with the y-axis and the slope of the fitted line, this class computes the:
			- squared pearson coefficient
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
			LinearRegression()
				: intercept_(0),
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

				@return If an error occured during the fit.

				@exception Exception::UnableToFit is thrown if fitting cannot be performed
			*/
    	template <typename Iterator>
			void computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin);

			/**
				@brief This function computes the best-fit linear regression coefficient \f$ (c_0) \f$
				of the model \f$ Y = c_1 X \f$ for the dataset \f$ (x, y) \f$.

				The values in x-dimension of the dataset \f$ (x,y) \f$ are given by the iterator range [x_begin,x_end)
				and the corresponding y-values start at position y_begin.

				For a  "x %" Confidence Interval use confidence_interval_P = x/100.
				For example the 95% Confidence Interval is supposed to be an interval that has a 95% chance of
				containing the true value of the parameter.

				@return If an error occured during the fit.

				@exception Exception::UnableToFit is thrown if fitting cannot be performed
			*/
    	template <typename Iterator>
			void computeRegressionNoIntercept(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin);

			/**
				@brief This function computes the best-fit linear regression coefficients \f$ (c_0,c_1) \f$
				of the model \f$ Y = c_0 + c_1 X \f$ for the weighted dataset \f$ (x, y) \f$.

				The values in x-dimension of the dataset \f$ (x, y) \f$ are given by the iterator range [x_begin,x_end)
				and the corresponding y-values start at position y_begin. They will be weighted by the
				values starting at w_begin.

				For a  "x %" Confidence Interval use confidence_interval_P = x/100.
				For example the 95% Confidence Interval is supposed to be an interval that has a 95% chance of
				containing the true value of the parameter.

				@return If an error occured during the fit.

				@exception Exception::UnableToFit is thrown if fitting cannot be performed
			*/
    	template <typename Iterator>
			void computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin);


			/// Non-mutable access to the y-intercept of the straight line
			DoubleReal getIntercept() const;
			/// Non-mutable access to the slope of the straight line
			DoubleReal getSlope() const;
			/// Non-mutable access to the x-intercept of the straight line
			DoubleReal getXIntercept() const;
			/// Non-mutable access to the lower border of confidence interval
			DoubleReal getLower() const;
			/// Non-mutable access to the upper border of confidence interval
			DoubleReal getUpper() const;
			/// Non-mutable access to the value of the t-distribution
			DoubleReal getTValue() const;
			/// Non-mutable access to the squared pearson coefficient
			DoubleReal getRSquared() const;
			/// Non-mutable access to the standard deviation of the residuals
			DoubleReal getStandDevRes() const;
			/// Non-mutable access to the residual mean
			DoubleReal getMeanRes() const;
			/// Non-mutable access to the standard error of the slope
			DoubleReal getStandErrSlope() const;
			/// Non-mutable access to the chi squared value
			DoubleReal getChiSquared() const;
			/// Non-mutable access to relelative standard deviation
			DoubleReal getRSD() const;

		 protected:

			/// The intercept of the fitted line with the y-axis
			double intercept_;
			/// The slope of the fitted line
			double slope_;
			/// The intercept of the fitted line with the x-axis
			double x_intercept_;
			/// The lower bound of the confidence intervall
			double lower_;
			/// The upper bound of the confidence intervall
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
			void computeGoodness_(double* X, double* Y, int N, double confidence_interval_P);

			/// Copies the distance(x_begin,x_end) elements starting at x_begin and y_begin into the arrays x_array and y_array
    	template <typename Iterator>
    	void iteratorRange2Arrays_(Iterator x_begin, Iterator x_end, Iterator y_begin, double* x_array,  double* y_array);

			/// Copy the distance(x_begin,x_end) elements starting at  x_begin, y_begin and w_begin into the arrays x_array, y_array and w_array
   	  template <typename Iterator>
			void iteratorRange3Arrays_ (Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin, double* x_array, double* y_array, double* w_array);

		private:

			/// Not implemented
			LinearRegression(const LinearRegression& arg);
			/// Not implemented
			LinearRegression& operator=(const LinearRegression& arg);
    };

    template <typename Iterator>
    void LinearRegression::computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin)
    {
      int N = int(distance(x_begin,x_end));

      double* X = new double[N];
      double* Y = new double[N];
      iteratorRange2Arrays_(x_begin,x_end,y_begin,X,Y);

      double cov00, cov01, cov11;

      // Compute the unweighted linear fit.
      // Get the intercept and the slope of the regression Y_hat=intercept_+slope_*X
      // and the value of Chi squared, the covariances of the intercept and the slope
      int error = gsl_fit_linear (X,1,Y,1,N,&intercept_, &slope_, &cov00, &cov01, &cov11, &chi_squared_);

      if (!error)
      {
        computeGoodness_(X,Y,N,confidence_interval_P);
      }

      delete [] X;
      delete [] Y;

			if (error)
			{
				throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"UnableToFit-LinearRegression","Could not fit a linear model to the data");
			}
    }

    template <typename Iterator>
    void LinearRegression::computeRegressionNoIntercept(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin)
    {
      int N = int(distance(x_begin,x_end));

      double* X = new double[N];
      double* Y = new double[N];
      iteratorRange2Arrays_(x_begin,x_end,y_begin,X,Y);

      double cov;

      // Compute the linear fit.
      // Get the intercept and the slope of the regression Y_hat=intercept_+slope_*X
      // and the value of Chi squared, the covariances of the intercept and the slope
      int error = gsl_fit_mul(X,1,Y,1,N,&slope_, &cov, &chi_squared_);
			intercept_ = 0.0;

      if (!error)
      {
        computeGoodness_(X,Y,N,confidence_interval_P);
      }

      delete [] X;
      delete [] Y;

			if (error)
			{
				throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"UnableToFit-LinearRegression","Could not fit a linear model to the data");
			}
    }

    template <typename Iterator>
    void LinearRegression::computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)
    {
      int N = int(distance(x_begin,x_end));

      double* X = new double[N];
      double* Y = new double[N];
      double* W = new double[N];
      iteratorRange3Arrays_(x_begin,x_end,y_begin,w_begin,X,Y,W);

      double cov00, cov01, cov11;

      // Compute the weighted linear fit.
      // Get the intercept and the slope of the regression Y_hat=intercept_+slope_*X
      // and the value of Chi squared, the covariances of the intercept and the slope
      int error = gsl_fit_wlinear (X,1,W,1,Y,1,N,&intercept_, &slope_, &cov00, &cov01, &cov11, &chi_squared_);

      if (!error)
      {
        computeGoodness_(X,Y,N,confidence_interval_P);
      }

      delete [] X;
      delete [] Y;
      delete [] W;

			if (error)
			{
				throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"UnableToFit-LinearRegression","Could not fit a linear model to the data");
			}
    }

    template <typename Iterator>
    void LinearRegression::iteratorRange2Arrays_(Iterator x_begin, Iterator x_end, Iterator y_begin, double* x_array, double* y_array)
    {
      int i=0;
      while (x_begin < x_end)
      {
        x_array[i]= *x_begin;
        y_array[i]= *y_begin;
        ++x_begin;
        ++y_begin;
        ++i;
      }
    }

    template <typename Iterator>
    void LinearRegression::iteratorRange3Arrays_(Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin, double* x_array, double* y_array, double* w_array)
    {
      int i=0;
      while (x_begin < x_end)
      {
        x_array[i]=*x_begin;
        y_array[i]=*y_begin;
        w_array[i]=*w_begin;
        ++x_begin;
        ++y_begin;
        ++w_begin;
        ++i;
      }
    }
  } // namespace Math
} // namespace OpenMS


#endif
