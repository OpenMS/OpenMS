// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

namespace OpenMS
{
	namespace Math
	{
		DoubleReal LinearRegression::getIntercept() const
		{
			return intercept_;
		}
		DoubleReal LinearRegression::getSlope() const
		{
			return slope_;
		}
		DoubleReal LinearRegression::getXIntercept() const
		{
			return x_intercept_;
		}
		DoubleReal LinearRegression::getLower() const
		{
			return lower_;
		}
		DoubleReal LinearRegression::getUpper() const
		{
			return upper_;
		}
		DoubleReal LinearRegression::getTValue() const
		{
			return t_star_;
		}
		DoubleReal LinearRegression::getRSquared() const
		{
			return r_squared_;
		}
		DoubleReal LinearRegression::getStandDevRes() const
		{
			return stand_dev_residuals_;
		}
			
		DoubleReal LinearRegression::getMeanRes() const
		{
			return mean_residuals_;
		}
			
		DoubleReal LinearRegression::getStandErrSlope() const
		{
			return stand_error_slope_;
		}
		DoubleReal LinearRegression::getChiSquared() const
		{
			return chi_squared_;
		}
			
		DoubleReal LinearRegression::getRSD() const
		{
			return rsd_;
		}

    void LinearRegression::computeGoodness_(double* X, double* Y,int N, double confidence_interval_P)
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

	}
}
