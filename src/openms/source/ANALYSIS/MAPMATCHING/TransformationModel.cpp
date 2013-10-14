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
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm>
#include <numeric>
#include <iterator>


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
                           ListUtils::create<String>("true,false"));
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
    x_.resize(mapping.size());
    y_.resize(mapping.size());
    size_t i = 0;
    for (map<DoubleReal, vector<DoubleReal> >::const_iterator it =
           mapping.begin(); it != mapping.end(); ++it, ++i)
    {
      x_[i] = it->first;
      // use average y value:
      y_[i] = accumulate(it->second.begin(), it->second.end(), 0.0) /
              it->second.size();
    }
    if (x_.size() < 3)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Cubic spline model needs at least 3 data points (with unique x values)");
    }
    interp_ = new Spline2d<double>(3, x_, y_);

    // linear model for extrapolation:
    TransformationModel::DataPoints lm_data(2);
    lm_data[0] = make_pair(x_.front(), y_.front());
    lm_data[1] = make_pair(x_.back(), y_.back());
    lm_ = new TransformationModelLinear(lm_data, Param());
  }

  TransformationModelInterpolated::~TransformationModelInterpolated()
  {
    delete interp_;
    delete lm_;
  }

  DoubleReal TransformationModelInterpolated::evaluate(const DoubleReal value)
  const
  {
    if ((value < x_.front()) || (value > x_.back()))     // extrapolate
    {
      return lm_->evaluate(value);
    }
    // interpolate:
    return interp_->eval(value);
  }

  void TransformationModelInterpolated::getDefaultParameters(Param & params)
  {
    params.clear();
    params.setValue("interpolation_type", "cspline",
                    "Type of interpolation to apply.");
    StringList types = ListUtils::create<String>("linear,polynomial,cspline,akima");
    params.setValidStrings("interpolation_type", types);
  }
}
