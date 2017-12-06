// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

// Spline2dInterpolator
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <map>
#include <vector>
#include <algorithm>
#include <numeric>

// AkimaInterpolator
#include <Wm5IntpAkimaNonuniform1.h>
#include <Wm5Math.h>

namespace OpenMS
{
  /**
   * @brief Spline2dInterpolator
   */
  class Spline2dInterpolator :
    public TransformationModelInterpolated::Interpolator
  {
public:
    Spline2dInterpolator() :
      spline_(nullptr)
    {
    }

    void init(std::vector<double>& x, std::vector<double>& y) override
    {
      // cleanup before we use a new one
      if (spline_ != (CubicSpline2d*) nullptr) delete spline_;

      // initialize spline
      spline_ = new CubicSpline2d(x, y);
    }

    double eval(const double& x) const override
    {
      return spline_->eval(x);
    }

    ~Spline2dInterpolator() override
    {
      delete spline_;
    }

private:
    CubicSpline2d* spline_;
    // Spline2d<double>* spline_;
  };

  /**
   * @brief AkimaInterpolator
   */
  class AkimaInterpolator :
    public TransformationModelInterpolated::Interpolator
  {
public:
    AkimaInterpolator() :
      interpolator_(nullptr)
    {}

    void init(std::vector<double>& x, std::vector<double>& y) override
    {
      if (interpolator_ != (Wm5::IntpAkimaNonuniform1<double>*) nullptr) delete interpolator_;
      // re-construct a new interpolator
      interpolator_ = new Wm5::IntpAkimaNonuniform1<double>(static_cast<int>(x.size()), &x.front(), &y.front());
    }

    double eval(const double& x) const override
    {
      return (* interpolator_)(x);
    }

    ~AkimaInterpolator() override
    {
      delete interpolator_;
    }

private:
    Wm5::IntpAkimaNonuniform1<double>* interpolator_;
  };

  /**
   * @brief LinearInterpolator.
   */
  class LinearInterpolator :
    public TransformationModelInterpolated::Interpolator
  {
public:
    LinearInterpolator()
    {}

    void init(std::vector<double>& x, std::vector<double>& y) override
    {
      // clear data
      x_.clear();
      y_.clear();

      // copy data
      // TODO: should we solve this using pointers to the original data?
      x_.insert(x_.begin(), x.begin(), x.end());
      y_.insert(y_.begin(), y.begin(), y.end());
    }

    double eval(const double& x) const override
    {
      // find nearest pair of points
      std::vector<double>::const_iterator it = std::upper_bound(x_.begin(), x_.end(), x);

      // interpolator is guaranteed to be only evaluated on points x, x_.front() =< x =< x x.back()
      // see TransformationModelInterpolated::evaluate

      // compute interpolation
      // the only point that is > then an element in our series is y_.back()
      // see call guarantee above
      if (it == x_.end())
      {
        return y_.back();
      }
      else
      {
        // interpolate .. invariant: idx > 0
        const SignedSize idx = it - x_.begin();
        const double x_0 = x_[idx - 1];
        const double x_1 = x_[idx];
        const double y_0 = y_[idx - 1];
        const double y_1 = y_[idx];

        return y_0 + (y_1 - y_0) * (x - x_0) / (x_1 - x_0);
      }
    }

    ~LinearInterpolator() override
    {
    }

private:
    /// x values
    std::vector<double> x_;
    /// y values
    std::vector<double> y_;
  };

  void TransformationModelInterpolated::preprocessDataPoints_(const DataPoints& data)
  {
    // need monotonically increasing x values (can't have the same value twice):
    std::map<double, std::vector<double> > mapping;
    for (TransformationModel::DataPoints::const_iterator it = data.begin();
         it != data.end();
         ++it)
    {
      mapping[it->first].push_back(it->second);
    }
    x_.resize(mapping.size());
    y_.resize(mapping.size());
    size_t i = 0;
    for (std::map<double, std::vector<double> >::const_iterator it = mapping.begin();
         it != mapping.end();
         ++it, ++i)
    {
      x_[i] = it->first;
      // use average y value:
      y_[i] = std::accumulate(it->second.begin(), it->second.end(), 0.0) / it->second.size();
    }

    // ensure that we have enough points for an interpolation
    if (x_.size() < 3)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Cubic spline model needs at least 3 data points (with unique x values)");
    }
  }

  TransformationModelInterpolated::TransformationModelInterpolated(const TransformationModel::DataPoints& data, const Param& params)
  {
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    // convert incoming data to x_ and y_
    preprocessDataPoints_(data);

    // choose the actual interpolation type
    const String interpolation_type = params_.getValue("interpolation_type");
    if (interpolation_type == "linear")
    {
      interp_ = new LinearInterpolator();
    }
    else if (interpolation_type == "cspline")
    {
      interp_ = new Spline2dInterpolator();
    }
    else if (interpolation_type == "akima")
    {
      interp_ = new AkimaInterpolator();
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "unknown/unsupported interpolation type '" + interpolation_type + "'");
    }

    // assign data
    interp_->init(x_, y_);

    // linear model for extrapolation:
    const String extrapolation_type = params_.getValue("extrapolation_type");
    if (extrapolation_type == "global-linear")
    {
      lm_front_ = new TransformationModelLinear(data, Param());
      lm_back_ = new TransformationModelLinear(data, Param());
    }
    else if (extrapolation_type == "two-point-linear")
    {
      TransformationModel::DataPoints lm_data(2);
      lm_data[0] = std::make_pair(x_.front(), y_.front());
      lm_data[1] = std::make_pair(x_.back(), y_.back()); // last point
      lm_front_ = new TransformationModelLinear(lm_data, Param());
      lm_back_ = new TransformationModelLinear(lm_data, Param());
    }
    else if (extrapolation_type == "four-point-linear")
    {
      TransformationModel::DataPoints lm_data(2);
      lm_data[0] = std::make_pair(x_[0], y_[0]); 
      lm_data[1] = std::make_pair(x_[1], y_[1]);
      lm_front_ = new TransformationModelLinear(lm_data, Param());

      lm_data[0] = std::make_pair(x_[ x_.size()-2 ], y_[ y_.size()-2] ); // second to last point
      lm_data[1] = std::make_pair(x_.back(), y_.back()); // last point
      lm_back_ = new TransformationModelLinear(lm_data, Param());
    }
    else
    {
      if (interp_) 
      {
        delete interp_;
      }

      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "unknown/unsupported extrapolation type '" + extrapolation_type + "'");
    }
  }

  TransformationModelInterpolated::~TransformationModelInterpolated()
  {
    if (interp_) delete interp_;
    if (lm_front_) delete lm_front_;
    if (lm_back_) delete lm_back_;
  }

  double TransformationModelInterpolated::evaluate(double value) const
  {
    if (value < x_.front()) // extrapolate front
    {
      return lm_front_->evaluate(value);
    }
    else if (value > x_.back()) // extrapolate back
    {
      return lm_back_->evaluate(value);
    }
    // interpolate:
    return interp_->eval(value);
  }

  void TransformationModelInterpolated::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("interpolation_type", "cspline", "Type of interpolation to apply.");
    StringList types = ListUtils::create<String>("linear,cspline,akima");
    params.setValidStrings("interpolation_type", types);
    params.setValue("extrapolation_type", "two-point-linear", "Type of extrapolation to apply: two-point-linear: use the first and last data point to build a single linear model, four-point-linear: build two linear models on both ends using the first two / last two points, global-linear: use all points to build a single linear model. Note that global-linear may not be continuous at the border.");
    StringList etypes = ListUtils::create<String>("two-point-linear,four-point-linear,global-linear");
    params.setValidStrings("extrapolation_type", etypes);
  }

} // namespace
