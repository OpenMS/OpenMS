// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_STATISTICS_LINEARREGRESSION_H
#define OPENMS_MATH_STATISTICS_LINEARREGRESSION_H

#include<OpenMS/CONCEPT/Types.h>

#include<iostream>
#include<vector>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_statistics.h>
#include<gsl/gsl_cdf.h>

namespace OpenMS
{
  namespace Math
  {
    /**
			@brief This class offers functions to perform least-squares fits to a straight line model, \f$ Y(c,x) = c_0 + c_1 x \f$.
	            
			It capsulates the GSL methods for a weighted and an unweighted linear regression.
	      
			Next to the intercept with the y-axis and the slope of the fitted line, this class
			computes the:
			<UL>
			<li> squared pearson coefficient</LI>       
			<li> value of the t-distribution</LI>      
			<li> standard deviation of the residuals</LI> 
			<li> standard error of the slope</LI>   
			<li> intercept with the x-axis (useful for additive series experiments)</LI>
			<li> lower border of confidence interval</LI>  
			<li> higher border of confidence interval</LI> 
			<li> chi squared value</LI>      
			<li> x mean</LI> 
			</UL>
			
			@todo Only the actual run-method should be a template (Clemens)
			
			@ingroup Math
    */
    template <typename Iterator>
    class LinearRegression
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
			{}

			/// Copy constructor.
			LinearRegression ( LinearRegression const & arg )
				: intercept_(arg.intercept_),
					slope_(arg.slope_),
					x_intercept_(arg.x_intercept_),
					lower_(arg.lower_),
					upper_(arg.upper_),
					t_star_(arg.t_star_),
					r_squared_(arg.r_squared_),
					stand_dev_residuals_(arg.stand_dev_residuals_),
					mean_residuals_(arg.mean_residuals_),
					stand_error_slope_(arg.stand_error_slope_),
					chi_squared_(arg.chi_squared_),
					rsd_(arg.rsd_)
			{}

			/// Assignment.
			LinearRegression & operator = ( LinearRegression const & arg )
			{
				// take care of self assignments
				if (this == &arg)
				{
					return *this;
				}

				intercept_ = arg.intercept_;
				slope_ = arg.slope_;
				x_intercept_ = arg.x_intercept_;
				lower_ = arg.lower_;
				upper_ = arg.upper_;
				t_star_ = arg.t_star_;
				r_squared_ = arg.r_squared_;
				stand_dev_residuals_ = arg.stand_dev_residuals_;
				mean_residuals_  = arg.mean_residuals_;
				stand_error_slope_ = arg.stand_error_slope_;
				chi_squared_ = arg.chi_squared_;
				rsd_ = arg.rsd_;

				return *this;
			}


			/// Destructor
			virtual ~LinearRegression()
			{}


			/**
			@brief This function computes the best-fit linear regression coefficients \f$ (c_0,c_1) \f$ 
			of the model \f$ Y = c_0 + c_1 X \f$ for the dataset \f$ (x, y) \f$.
                    
			The values in x-dimension of the dataset \f$ (x,y) \f$ are given by the iterator range [x_begin,x_end)
			and the corresponding y-values start at position y_begin.
                  

			For a  "x %" Confidence Interval use confidence_interval_P = x/100.
			For example the 95% Confidence Interval is supposed to be an interval that has a 95% chance of 
			containing the true value of the parameter.
			*/
			int computeRegression
			(double confidence_interval_P,
			 Iterator x_begin,
			 Iterator x_end,
			 Iterator y_begin);

			/**
			@brief This function computes the best-fit linear regression coefficients \f$ (c_0,c_1) \f$ 
			of the model \f$ Y = c_0 + c_1 X \f$ for the weighted dataset \f$ (x, y) \f$.
                    
			The values in x-dimension of the dataset \f$ (x, y) \f$ are given by the iterator range [x_begin,x_end)
			and the corresponding y-values start at position y_begin. They will be weighted by the 
			values starting at w_begin.

			For a  "x %" Confidence Interval use confidence_interval_P = x/100.
			For example the 95% Confidence Interval is supposed to be an interval that has a 95% chance of 
			containing the true value of the parameter.
			*/
			int computeRegressionWeighted
			(double confidence_interval_P,
			 Iterator x_begin,
			 Iterator x_end,
			 Iterator y_begin,
			 Iterator w_begin);


			/// Non-mutable access to the y-intercept of the straight line
			inline DoubleReal getIntercept() const
			{
				return intercept_;
			}
			/// Non-mutable access to the slope of the straight line
			inline DoubleReal getSlope() const
			{
				return slope_;
			}
			/// Non-mutable access to the x-intercept of the straight line
			inline DoubleReal getXIntercept() const
			{
				return x_intercept_;
			}
			/// Non-mutable access to the lower border of confidence interval
			inline DoubleReal getLower() const
			{
				return lower_;
			}
			/// Non-mutable access to the upper border of confidence interval
			inline DoubleReal getUpper() const
			{
				return upper_;
			}
			/// Non-mutable access to the value of the t-distribution
			inline DoubleReal getTValue() const
			{
				return t_star_;
			}
			/// Non-mutable access to the squared pearson coefficient
			inline DoubleReal getRSquared() const
			{
				return r_squared_;
			}
			/// Non-mutable access to the standard deviation of the residuals
			inline DoubleReal getStandDevRes() const
			{
				return stand_dev_residuals_;
			}
				
			/// Non-mutable access to the residual mean 
			inline DoubleReal getMeanRes() const
			{
				return mean_residuals_;
			}
				
			/// Non-mutable access to the standard error of the slope
			inline DoubleReal getStandErrSlope() const
			{
				return stand_error_slope_;
			}
			/// Non-mutable access to the chi squared value
			inline DoubleReal getChiSquared() const
			{
				return chi_squared_;
			}
				
			/// Non-mutable access to relelative standard deviation
			inline DoubleReal getRSD() const
			{
				return rsd_;
			}


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
			void computeGoodness_
			(double* X,
			 double* Y,
			 int N,
			 double confidence_interval_P);

			/// Copies the distance(x_begin,x_end) elements starting at x_begin and y_begin into the arrays x_array and y_array
			void iteratorRange2Arrays_
			(Iterator x_begin,
			 Iterator x_end,
			 Iterator y_begin,
			 double* x_array,
			 double* y_array);

			/** @brief Copy the distance(x_begin,x_end) elements starting at
			x_begin, y_begin and w_begin into the arrays x_array, 
			y_array and w_array
			*/
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
    int LinearRegression<Iterator>::computeRegression
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
    int LinearRegression<Iterator>::computeRegressionWeighted
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

      // Mean of abscissa and ordinate values
      double x_mean = gsl_stats_mean(X,1,N);
			double y_mean = gsl_stats_mean(Y,1,N);
				
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
      double sum = 0;
      for (Int i = 0; i < N; ++i)
      {
        double x_i = fabs(Y[i] - (intercept_ + slope_ * X[i]));
        sum += x_i;
      }
			mean_residuals_       = sum / N;
      stand_dev_residuals_ = sqrt((chi_squared_ - (sum*sum)/N)/(N-1));

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
			
			double tmp = 0;
			for (Int i = 0; i<N;++i)
			{
				tmp += (X[i] - x_mean) * (X[i] - x_mean);
			}
			
			// 			cout << "100.0 / abs( x_intercept_ ) " << (100.0 / fabs( x_intercept_ )) << endl;
			// 			cout << "tmp : " << tmp << endl;
			// 			cout << "slope_ " << slope_ << endl;
			// 			cout << "y_mean " << y_mean << endl;
			// 			cout << "N " << N << endl;
			// 			cout << "stand_dev_residuals_ " << stand_dev_residuals_ << endl;
			// 			cout << " (1.0/ (double) N)  " <<  (1.0/ (double) N)  << endl;
			// 			cout << "sx hat " << (stand_dev_residuals_ / slope_) * sqrt(  (1.0/ (double) N) * (y_mean / (slope_ * slope_ * tmp ) ) ) << endl;
	
			// compute relative standard deviation (non-standard formula, taken from Mayr et al. (2006) )
			rsd_ = (100.0 / fabs( x_intercept_ )) * (stand_dev_residuals_ / slope_) * sqrt(  (1.0/ (double) N) * (y_mean / (slope_ * slope_ * tmp ) ) ) ; 
			
			if (rsd_ < 0.0)
			{
				std::cout << "rsd < 0.0 " << std::endl;
				std::cout <<   "Intercept                                " << intercept_
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
									<< std::endl; 	
			}
			
 
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
  } // namespace Math
} // namespace OpenMS


#endif
