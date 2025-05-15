// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

// Spline2dInterpolator
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <numeric>

// AkimaInterpolator
#include <Mathematics/IntpAkimaNonuniform1.h>

namespace OpenMS
{
  /**
   * @brief Spline2dInterpolator
   */
  class Spline2dInterpolator :
    public TransformationModelInterpolated::Interpolator
  {
public:
    Spline2dInterpolator() 
      
    = default;

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
    CubicSpline2d* spline_{nullptr};
    // Spline2d<double>* spline_;
  };

  /**
   * @brief AkimaInterpolator
   */
  class AkimaInterpolator :
    public TransformationModelInterpolated::Interpolator
  {
public:
    AkimaInterpolator() 
      
    = default;

    void init(std::vector<double>& x, std::vector<double>& y) override
    {
      if (interpolator_ != (gte::IntpAkimaNonuniform1<double>*) nullptr) delete interpolator_;
      // re-construct a new interpolator
      interpolator_ = new gte::IntpAkimaNonuniform1<double>(static_cast<int>(x.size()), &x.front(), &y.front());
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
    gte::IntpAkimaNonuniform1<double>* interpolator_{nullptr};
  };

  /**
   * @brief LinearInterpolator.
   */
  class LinearInterpolator :
    public TransformationModelInterpolated::Interpolator
  {
public:
    LinearInterpolator()
    = default;

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
    = default;

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

  void TransformationModelInterpolated::preprocessDataPoints_(const std::vector<std::pair<double,double>>& data)
  {
    // need monotonically increasing x values (can't have the same value twice):
    std::map<double, std::vector<double> > mapping;
    for (std::vector<std::pair<double,double>>::const_iterator it = data.begin();
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

  TransformationModelInterpolated::TransformationModelInterpolated(const std::vector<std::pair<double,double>>& data, const Param& params, bool preprocess = true)
  {
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    // convert incoming data to x_ and y_
    if (preprocess)
    {
      preprocessDataPoints_(data);
    }
    else
    {
      x_.resize(data.size());
      y_.resize(data.size());
      for (const std::pair<double,double>& pair : data)
      {
        x_.push_back(pair.first);
        y_.push_back(pair.second);
      }
    }


    // choose the actual interpolation type
    const String interpolation_type = params_.getValue("interpolation_type").toString();
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
    const String extrapolation_type = params_.getValue("extrapolation_type").toString();
    if (extrapolation_type == "global-linear")
    {
      std::vector<TransformationModel::DataPoint> bloated_data{};
      bloated_data.resize(x_.size());
      //uff... well here we go.. adding an empty string
      for (Size s = 0; s < x_.size(); ++s)
      {
        bloated_data.emplace_back(TransformationModel::DataPoint(x_[s],y_[s]));
      }
      lm_front_ = new TransformationModelLinear(bloated_data, Param());
      lm_back_ = new TransformationModelLinear(bloated_data, Param());
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

  TransformationModelInterpolated::TransformationModelInterpolated(const TransformationModel::DataPoints& data, const Param& params)
  {
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    // convert incoming data to x_ and y_
    preprocessDataPoints_(data);

    // choose the actual interpolation type
    const String interpolation_type = params_.getValue("interpolation_type").toString();
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
    const String extrapolation_type = params_.getValue("extrapolation_type").toString();
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
    params.setValidStrings("interpolation_type", {"linear","cspline","akima"});
    params.setValue("extrapolation_type", "two-point-linear", "Type of extrapolation to apply: two-point-linear: use the first and last data point to build a single linear model, four-point-linear: build two linear models on both ends using the first two / last two points, global-linear: use all points to build a single linear model. Note that global-linear may not be continuous at the border.");
    params.setValidStrings("extrapolation_type", {"two-point-linear","four-point-linear","global-linear"});
  }

} // namespace
