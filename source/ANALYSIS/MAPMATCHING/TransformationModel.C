// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <numeric>

using namespace std;

namespace OpenMS
{
  TransformationModelLinear::TransformationModelLinear(
    const TransformationModel::DataPoints & data, const Param & params)
  {
    params_ = params;
    data_given_ = !data.empty();
    if (!data_given_ && params.exists("slope") && (params.exists("intercept")))
    {
      // don't estimate parameters, use given values
      slope_ = params.getValue("slope");
      intercept_ = params.getValue("intercept");
    }
    else     // estimate parameters from data
    {
      Param defaults;
      getDefaultParameters(defaults);
      params_.setDefaults(defaults);
      symmetric_ = params_.getValue("symmetric_regression") == "true";

      size_t size = data.size();
      if (size == 0)       // no data
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                         "no data points for 'linear' model");
      }
      else if (size == 1)       // degenerate case, but we can still do something
      {
        slope_ = 1.0;
        intercept_ = data[0].second - data[0].first;
      }
      else       // compute least-squares fit
      {
        vector<double> x(size), y(size);
        for (size_t i = 0; i < size; ++i)
        {
          if (symmetric_)
          {
            x[i] = data[i].second + data[i].first;
            y[i] = data[i].second - data[i].first;
          }
          else
          {
            x[i] = data[i].first;
            y[i] = data[i].second;
          }
        }
        double cov00, cov01, cov11, sumsq;         // covariance values, sum of squares
        double * x_start = &(x[0]), * y_start = &(y[0]);
        deprecated_gsl_fit_linear(x_start, 1, y_start, 1, size, &intercept_, &slope_,
                       &cov00, &cov01, &cov11, &sumsq);

        if (symmetric_)         // undo coordinate transformation:
        {
          slope_ = (1.0 + slope_) / (1.0 - slope_);
          intercept_ = intercept_ * 1.41421356237309504880;           // 1.41... = sqrt(2)
        }
      }
    }
  }

  TransformationModelLinear::~TransformationModelLinear()
  {
  }

  DoubleReal TransformationModelLinear::evaluate(const DoubleReal value) const
  {
    return slope_ * value + intercept_;
  }

  void TransformationModelLinear::invert()
  {
    if (slope_ == 0)
      throw Exception::DivisionByZero(__FILE__, __LINE__,
                                      __PRETTY_FUNCTION__);
    intercept_ = -intercept_ / slope_;
    slope_ = 1.0 / slope_;
    // update parameters:
    if (params_.exists("slope") && (params_.exists("intercept")))
    {
      params_.setValue("slope", slope_, params_.getDescription("slope"));
      params_.setValue("intercept", intercept_,
                       params_.getDescription("intercept"));
    }
  }

  void TransformationModelLinear::getParameters(DoubleReal & slope,
                                                DoubleReal & intercept) const
  {
    slope = slope_;
    intercept = intercept_;
  }

  void TransformationModelLinear::getDefaultParameters(Param & params)
  {
    params.clear();
    params.setValue("symmetric_regression", "false", "Perform linear regression"
                                                     " on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.");
    params.setValidStrings("symmetric_regression",
                           StringList::create("true,false"));
  }

  TransformationModelInterpolated::TransformationModelInterpolated(
    const TransformationModel::DataPoints & data, const Param & params)
  {
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    // need monotonically increasing x values (can't have the same value twice):
    map<DoubleReal, vector<DoubleReal> > mapping;
    for (TransformationModel::DataPoints::const_iterator it = data.begin();
         it != data.end(); ++it)
    {
      mapping[it->first].push_back(it->second);
    }
    size_ = mapping.size();
    x_.resize(size_);
    y_.resize(size_);
    size_t i = 0;
    for (map<DoubleReal, vector<DoubleReal> >::const_iterator it =
           mapping.begin(); it != mapping.end(); ++it, ++i)
    {
      x_[i] = it->first;
      // use average y value:
      y_[i] = accumulate(it->second.begin(), it->second.end(), 0.0) /
              it->second.size();
    }

    String interpolation_type = params_.getValue("interpolation_type");
    const deprecated_gsl_interp_type * type;
    if (interpolation_type == "linear")
      type = deprecated_wrapper_get_gsl_interp_linear();
    else if (interpolation_type == "polynomial")
      type = deprecated_wrapper_get_gsl_interp_polynomial();
    else if (interpolation_type == "cspline")
      type = deprecated_wrapper_get_gsl_interp_cspline();
    else if (interpolation_type == "akima")
      type = deprecated_wrapper_get_gsl_interp_akima();
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "unknown/unsupported interpolation type '" + interpolation_type + "'");
    }

    size_t min_size = deprecated_wrapper_gsl_interp_type_get_min_size( type );
    if (size_ < min_size)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'" + interpolation_type + "' interpolation model needs at least " + String(min_size) + " data points (with unique x values)");
    }

    interp_ = deprecated_gsl_interp_alloc(type, size_);
    acc_ = deprecated_gsl_interp_accel_alloc();
    double * x_start = &(x_[0]), * y_start = &(y_[0]);
    deprecated_gsl_interp_init(interp_, x_start, y_start, size_);

    // linear model for extrapolation:
    TransformationModel::DataPoints lm_data(2);
    lm_data[0] = make_pair(x_[0], y_[0]);
    lm_data[1] = make_pair(x_[size_ - 1], y_[size_ - 1]);
    lm_ = new TransformationModelLinear(lm_data, Param());
  }

  TransformationModelInterpolated::~TransformationModelInterpolated()
  {
    deprecated_gsl_interp_free(interp_);
    deprecated_gsl_interp_accel_free(acc_);
    delete lm_;
  }

  DoubleReal TransformationModelInterpolated::evaluate(const DoubleReal value)
  const
  {
    if ((value < x_[0]) || (value > x_[size_ - 1]))     // extrapolate
    {
      return lm_->evaluate(value);
    }
    // interpolate:
    const double * x_start = &(x_[0]), * y_start = &(y_[0]);
    return deprecated_gsl_interp_eval(interp_, x_start, y_start, value, acc_);
  }

  void TransformationModelInterpolated::getDefaultParameters(Param & params)
  {
    params.clear();
    params.setValue("interpolation_type", "cspline",
                    "Type of interpolation to apply.");
    StringList types = StringList::create("linear,polynomial,cspline,akima");
    params.setValidStrings("interpolation_type", types);
  }

  TransformationModelBSpline::TransformationModelBSpline(
    const TransformationModel::DataPoints & data, const Param & params)
  {
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    if (data.size() < 4)     // TODO: check number
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'b_spline' model needs at least four data points");
    }
    Size num_breakpoints = params_.getValue("num_breakpoints");
    String break_positions = params_.getValue("break_positions");
    if ((break_positions != "uniform") && (break_positions != "quantiles"))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "parameter 'break_positions' for 'b_spline' model must be 'uniform' or 'quantiles'");
    }

    size_ = data.size();
    x_ = deprecated_gsl_vector_alloc(size_);
    y_ = deprecated_gsl_vector_alloc(size_);
    w_ = deprecated_gsl_vector_alloc(size_);
    for (size_t i = 0; i < size_; ++i)
    {
      deprecated_gsl_vector_set(x_, i, data[i].first);
      deprecated_gsl_vector_set(y_, i, data[i].second);
      deprecated_gsl_vector_set(w_, i, 1.0);       // TODO: non-uniform weights
    }
    deprecated_gsl_vector_minmax(x_, &xmin_, &xmax_);

    // set up cubic (k = 4) spline workspace:
    if (num_breakpoints < 2)
    {
      num_breakpoints = 2;
      LOG_WARN << "Warning: Increased parameter 'num_breakpoints' to 2 (minimum)." << endl;
    }
    else if (num_breakpoints > size_ - 2)
    {
      num_breakpoints = size_ - 2;
      LOG_WARN << "Warning: Decreased parameter 'num_breakpoints' to " + String(num_breakpoints) + " (maximum for this number of data points)." << endl;
    }
    workspace_ = deprecated_gsl_bspline_alloc(4, num_breakpoints);
    if (break_positions == "uniform")
    {
      deprecated_gsl_bspline_knots_uniform(xmin_, xmax_, workspace_);
    }
    else
    {
      vector<double> quantiles(num_breakpoints, 1.0);
      double step = 1.0 / (num_breakpoints - 1);
      for (Size i = 0; i < num_breakpoints - 1; ++i)
      {
        quantiles[i] = i * step;
      }
      deprecated_gsl_vector * breakpoints;
      breakpoints = deprecated_gsl_vector_alloc(num_breakpoints);
      getQuantiles_(x_, quantiles, breakpoints);
      deprecated_gsl_bspline_knots(breakpoints, workspace_);
      deprecated_gsl_vector_free(breakpoints);
    }
    ncoeffs_ = deprecated_gsl_bspline_ncoeffs(workspace_);
    deprecated_gsl_vector_minmax(
    		deprecated_wrapper_gsl_bspline_workspace_get_knots(workspace_), &xmin_, &xmax_);
    computeFit_();
  }

  TransformationModelBSpline::~TransformationModelBSpline()
  {
    deprecated_gsl_bspline_free(workspace_);
    deprecated_gsl_vector_free(bsplines_);
    deprecated_gsl_vector_free(coeffs_);
    deprecated_gsl_matrix_free(cov_);
    deprecated_gsl_vector_free(x_);
    deprecated_gsl_vector_free(y_);
    deprecated_gsl_vector_free(w_);
  }

  void TransformationModelBSpline::getQuantiles_(
    const deprecated_gsl_vector * x, const vector<double> & quantiles, deprecated_gsl_vector * results)
  {
    deprecated_gsl_vector * x_sort;
    x_sort = deprecated_gsl_vector_alloc(
    		deprecated_wrapper_gsl_vector_get_size(x) );
    deprecated_gsl_vector_memcpy(x_sort, x);
    deprecated_gsl_sort_vector(x_sort);
    for (Size i = 0; i < quantiles.size(); ++i)
    {
      double q = deprecated_gsl_stats_quantile_from_sorted_data(
    		  deprecated_wrapper_gsl_vector_get_data(x_sort), 1,
    		  deprecated_wrapper_gsl_vector_get_size(x),quantiles[i]);
      deprecated_gsl_vector_set(results, i, q);
    }
    deprecated_gsl_vector_free(x_sort);
  }

  void TransformationModelBSpline::computeFit_()
  {
    // construct the fit matrix:
    deprecated_gsl_matrix * fit_matrix = deprecated_gsl_matrix_alloc(size_, ncoeffs_);
    bsplines_ = deprecated_gsl_vector_alloc(ncoeffs_);
    for (size_t i = 0; i < size_; ++i)
    {
      double xi = deprecated_gsl_vector_get(x_, i);
      deprecated_gsl_bspline_eval(xi, bsplines_, workspace_);
      for (size_t j = 0; j < ncoeffs_; ++j)
      {
        double bspline = deprecated_gsl_vector_get(bsplines_, j);
        deprecated_gsl_matrix_set(fit_matrix, i, j, bspline);
      }
    }
    // do the fit:
    deprecated_gsl_multifit_linear_workspace * multifit = deprecated_gsl_multifit_linear_alloc(
      size_, ncoeffs_);
    coeffs_ = deprecated_gsl_vector_alloc(ncoeffs_);
    cov_ = deprecated_gsl_matrix_alloc(ncoeffs_, ncoeffs_);
    double chisq;
    deprecated_gsl_multifit_wlinear(fit_matrix, w_, y_, coeffs_, cov_, &chisq, multifit);
    // clean-up:
    deprecated_gsl_matrix_free(fit_matrix);
    deprecated_gsl_multifit_linear_free(multifit);
    // for linear extrapolation (natural spline):
    computeLinear_(xmin_, slope_min_, offset_min_, sd_err_left_);
    computeLinear_(xmax_, slope_max_, offset_max_, sd_err_right_);
  }

  void TransformationModelBSpline::computeLinear_(
    const double pos, double & slope, double & offset, double & sd_err)
  {
    deprecated_gsl_bspline_deriv_workspace * deriv_workspace = deprecated_gsl_bspline_deriv_alloc(4);
    deprecated_gsl_matrix * deriv = deprecated_gsl_matrix_alloc(ncoeffs_, 2);
    deprecated_gsl_bspline_deriv_eval(pos, 1, deriv, workspace_, deriv_workspace);
    deprecated_gsl_bspline_deriv_free(deriv_workspace);
    double results[2];
    for (size_t j = 0; j < 2; ++j)
    {
      for (size_t i = 0; i < ncoeffs_; ++i)
      {
        deprecated_gsl_vector_set(bsplines_, i, deprecated_gsl_matrix_get(deriv, i, j));
      }
      deprecated_gsl_multifit_linear_est(bsplines_, coeffs_, cov_, &results[j], &sd_err);
      // TODO: find a good way to estimate the error
      // e.g. by reducing "num_breakpoints" and comparing to current?!
      if (sd_err >= 1)
      {
        // this SD seems on a much smaller scale than "yerr" in "evaluate"...
        // see GSL's "multilinear.c::gsl_multifit_linear_est" for details
        LOG_WARN << "Warning: B-spline extrapolation is unreliable. Consider reducing 'num_breakpoints' or using another model." << endl;
      }
    }
    deprecated_gsl_matrix_free(deriv);
    offset = results[0];
    slope = results[1];
  }

  DoubleReal TransformationModelBSpline::evaluate(const DoubleReal value) const
  {
    DoubleReal result;
    if (value < xmin_)     // extrapolate on left side
    {
      result = offset_min_ - slope_min_ * (xmin_ - value);
    }
    else if (value > xmax_)     // extrapolate on right side
    {
      result = offset_max_ + slope_max_ * (value - xmax_);
    }
    else     // evaluate B-splines
    {
      double yerr;
      deprecated_gsl_bspline_eval(value, bsplines_, workspace_);
      deprecated_gsl_multifit_linear_est(bsplines_, coeffs_, cov_, &result, &yerr);
    }
    return result;
  }

  void TransformationModelBSpline::getDefaultParameters(Param & params)
  {
    params.clear();
    params.setValue("num_breakpoints", 5, "Number of breakpoints of the cubic spline in the smoothing step. More breakpoints mean less smoothing. Reduce this number if the transformation has an unexpected shape.");
    params.setMinInt("num_breakpoints", 2);
    params.setValue("break_positions", "uniform", "How to distribute the breakpoints on the retention time scale. 'uniform': intervals of equal size; 'quantiles': equal number of data points per interval.");
    params.setValidStrings("break_positions", StringList::create("uniform,quantiles"));
  }

}
