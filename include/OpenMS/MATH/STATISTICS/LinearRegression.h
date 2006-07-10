// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: LinearRegression.h,v 1.10 2006/03/28 16:19:59 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_MATH_STATISTICS_LINEARREGRESSION_H
#define OPENMS_MATH_STATISTICS_LINEARREGRESSION_H

#include<iostream>
#include<vector> 
#include<gsl/gsl_fit.h>
#include<gsl/gsl_statistics.h>
#include<gsl/gsl_cdf.h>

namespace OpenMS
{
	/**
		@brief Linear regression 
		
		@todo document me, add to Math namespace (Eva)
		
		@ingroup Math
	*/
	template <typename Iterator> 
	class LinearRegression
  {
	 public:

    LinearRegression(){}
    virtual ~LinearRegression(){}
    
	   
		/** @name Computing the intercept of a linear regression with the x-axis
		 */
		//@{	
		// The 95% Confidence Interval is supposed to be an interval that has a 95% chance of 
		// containing the true value of the parameter.
    // To Compute the "x %" Confidence Interval use confidence_interval_P= x/100
		/// Compute the x-intercept of the weighted linear regression
    int computeInterceptXAxisWeighted
		(double confidence_interval_P,
		 Iterator x_begin,
		 Iterator x_end,
		 Iterator y_begin,
		 Iterator w_begin);
		
		/// Compute the x-intercept of the unweighted linear regression	
    int computeInterceptXAxis
		(double confidence_interval_P,
		 Iterator x_begin,
		 Iterator x_end,
		 Iterator y_begin);

		/**	@name Accessors */
		//@{
		inline const double& getIntercept() const { return intercept_; }
		///
		inline const double& getSlope() const { return slope_; }
		///
		inline const double& getXIntercept() const { return x_intercept_; }
		///
		inline const double& getLower() const { return lower_; }
		///
		inline const double& getUpper() const { return upper_; }
		///
		inline const double& getTValue() const { return t_star_; }
		///
		inline const double& getRSquared() const { return r_squared_; }
		///
		inline const double& getStandDevRes() const { return stand_dev_residuals_; }
		///
		inline const double& getStandErrSlope() const { return stand_error_slope_; }
		///
		inline const double& getChiSquared() const { return chi_squared_; }
		//@}
		
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
		/// The value of the t-distribution
    double t_star_;
		/// The squared correlation coefficient (Pearson)
    double r_squared_;
		/// The standard deviation of the ressiduals
    double stand_dev_residuals_;
		/// The standard error of the slope
    double stand_error_slope_;
    /// The value of the Chi Squared 
    double chi_squared_;  


    // Compute the goodness of the fitted regression line
    void computeGoodness_
		(double* X, 
		 double* Y,
		 int N, 
		 double confidence_interval_P);
	
		// Copy the distance(x_begin,x_end) elements starting at
		// x_begin and y_begin into the arrays x_array and y_array
    void iteratorRange2Arrays_
		(Iterator x_begin, 
		 Iterator x_end, 
		 Iterator y_begin, 
		 double* x_array, 
		 double* y_array);
     
		// Copy the distance(x_begin,x_end) elements starting at
		// x_begin, y_begin and w_begin into the arrays x_array, 
		// y_array and w_array
    void iteratorRange3Arrays_
		(Iterator x_begin, 
		 Iterator x_end, 
		 Iterator y_begin, 
		 Iterator w_begin, 
		 double* x_array, 
		 double* y_array, 
		 double* w_array);
   
	};


	template <typename Iterator> 
	int LinearRegression<Iterator>::computeInterceptXAxis
	(double confidence_interval_P,
	 Iterator x_begin,
	 Iterator x_end,
	 Iterator y_begin
	)
	{
		int N=distance(x_begin,x_end);
				
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
		
		return error;
	}


	template <typename Iterator> 
	int LinearRegression<Iterator>::computeInterceptXAxisWeighted
	(double confidence_interval_P,
	 Iterator x_begin,
	 Iterator x_end,
	 Iterator y_begin,
	 Iterator w_begin
	)
	{
		int N=distance(x_begin,x_end);
				
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
		
		return error;
	}

	 
	template <typename Iterator> 
	void LinearRegression<Iterator>::computeGoodness_(double* X, double* Y,int N, double confidence_interval_P)
	{
		// Variance and Covariances
		double var_X = gsl_stats_variance(X,1,N);
		double var_Y = gsl_stats_variance(Y,1,N);
		double cov_XY = gsl_stats_covariance(X,1,Y,1,N);
			
		// Mean of X
		double x_mean = gsl_stats_mean(X,1,N);
			
		// S_xx
		double s_XX = 0;
		for (int i=0; i<N; ++i)
		{
			double d=(X[i]-x_mean);
			s_XX +=d*d;
		}
			
		// Compute the squared Pearson coefficient
		r_squared_=(cov_XY*cov_XY)/(var_X*var_Y);
			
		// The standard deviation of the residuals
		stand_dev_residuals_=sqrt(chi_squared_/(N-2));
			
		// The Standard error of the slope
		stand_error_slope_=stand_dev_residuals_/sqrt(s_XX);
			
		// and the intersection of Y_hat with the x-axis
		x_intercept_= -(intercept_/slope_);	
			
		double P=1-(1-confidence_interval_P)/2;
		t_star_=gsl_cdf_tdist_Pinv(P,N-2);
			
		//Compute the asymmetric 95% confidence intervall of around the X-intercept
		double g =(t_star_/(slope_/stand_error_slope_));
		g *= g;
		double left = (x_intercept_-x_mean)*g;
		double bottom = 1-g;
		double d=(x_intercept_-x_mean);
		double right = t_star_*(stand_dev_residuals_/slope_)*sqrt((d*d)/s_XX + (bottom/N));
						
		// Confidence intervall lower_ <= X_intercept <= upper_ 
		lower_ = x_intercept_+(left+right)/bottom;
		upper_ = x_intercept_+(left-right)/bottom;
			
		if (lower_ > upper_)
		{
			std::swap(lower_,upper_);
		}

		// 
		/*std::cerr <<   "Intercept                                " << intercept_
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
							<< "\nstand_error_slope/slope_                 " << (stand_dev_residuals_/slope_)
							<< "\nCoefficient of Variation                 " << (stand_dev_residuals_/slope_)/x_mean*100  << std::endl
							<< "========================================="
							<< std::endl;	*/			
	}
		


	template <typename Iterator> 
	void LinearRegression<Iterator>::iteratorRange2Arrays_
	(Iterator x_begin, 
	 Iterator x_end, 
	 Iterator y_begin, 
	 double* x_array, 
	 double* y_array)
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
	void LinearRegression<Iterator>::iteratorRange3Arrays_
	(Iterator x_begin, 
	 Iterator x_end, 
	 Iterator y_begin, 
	 Iterator w_begin, 
	 double* x_array, 
	 double* y_array, 
	 double* w_array)
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
}


#endif
