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
// $Maintainer: Timo Sachsenberg  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/distributions.hpp>

using boost::math::detail::inverse_students_t;

namespace OpenMS
{
  namespace Math
  {
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

    void LinearRegression::computeGoodness_(const std::vector<Wm5::Vector2d>& points, double confidence_interval_P)
    {
      OPENMS_PRECONDITION(static_cast<unsigned>(points.size()) > 2, 
          "Cannot compute goodness of fit for regression with less than 3 data points");
      // specifically, boost throws an exception for a t-distribution with zero df

      unsigned N = static_cast<unsigned>(points.size());
      std::vector<double> X; X.reserve(N);
      std::vector<double> Y; Y.reserve(N);
      for (unsigned i = 0; i < N; ++i)
      {
        X.push_back(points[i].X());
        Y.push_back(points[i].Y());
      }

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
                  << "\nstand_error_slope/slope_                 " << (stand_dev_residuals_ / slope_)
                  << "\nCoefficient of Variation                 " << (stand_dev_residuals_ / slope_) / x_mean * 100  << std::endl
                  << "========================================="
                  << std::endl;
      }
    }

  }
}

