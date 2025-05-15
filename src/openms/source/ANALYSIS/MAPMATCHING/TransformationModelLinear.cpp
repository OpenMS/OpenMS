// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <Mathematics/Vector2.h>
#include <Mathematics/ApprHeightLine2.h>

namespace OpenMS
{

  TransformationModelLinear::TransformationModelLinear(const TransformationModel::DataPoints& data, const Param& params) :
    TransformationModel(data, params) // initializes model
  {
    data_given_ = !data.empty();

    if (!data_given_ && params.exists("slope") && params.exists("intercept"))
    {
      // don't estimate parameters, use given values
      slope_ = params.getValue("slope");
      intercept_ = params.getValue("intercept");
    }
    else // estimate parameters from data
    {
      Param defaults;
      getDefaultParameters(defaults);
      params_.setDefaults(defaults);
      symmetric_ = params_.getValue("symmetric_regression") == "true";
      // weight the data (if weighting is specified)
      TransformationModel::DataPoints data_weighted = data;
      // TrafoXML's prior to OpenMS 3.0 have x/y_weight = "" if unweighted 
      if ((params.exists("x_weight") && params.getValue("x_weight") != "x" && params.getValue("x_weight") != "") ||
          (params.exists("y_weight") && params.getValue("y_weight") != "y" && params.getValue("y_weight") != ""))
      {
        weightData(data_weighted);
      }

      size_t size = data_weighted.size();
      std::vector<gte::Vector2<double>> points;
      if (size == 0) // no data
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "no data points for 'linear' model");
      }
      else if (size == 1) // degenerate case, but we can still do something
      {
        slope_ = 1.0;
        intercept_ = data_weighted[0].second - data_weighted[0].first;
      }
      else if (size == 2)
      {
        // if the two points are too close, gte::HeightLineFit2 can't fit a line
        // but in the special case of two points, there is an exact solution and we don't need a least-sqaures fit
        slope_ = (data_weighted[1].second - data_weighted[0].second) / (data_weighted[1].first - data_weighted[0].first);
        intercept_ = data_weighted[0].second - (slope_ * data_weighted[0].first);

        if (std::isinf(slope_) || std::isnan(slope_) || std::isinf(intercept_) || std::isnan(intercept_))
        {
          throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TransformationModelLinear", "Unable to fit linear transformation to the two data points.");
        }
      }
      else // compute least-squares fit
      {
        for (size_t i = 0; i < size; ++i)
        {
          points.emplace_back(std::initializer_list<double>{data_weighted[i].first, data_weighted[i].second});
        }
        auto line = gte::ApprHeightLine2<double>();
        if (!line.Fit(static_cast<int>(size), &points.front()))
        {
          throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TransformationModelLinear", "Unable to fit linear transformation to data points.");
        }
      slope_ = line.GetParameters().second[0];
      intercept_ = -slope_ * line.GetParameters().first[0] + line.GetParameters().first[1];
      }
      // update params
      params_.setValue("slope", slope_);
      params_.setValue("intercept", intercept_);
    }
  }

  double TransformationModelLinear::evaluate(double value) const
  {
    if (!weighting_)
    {
      return slope_ * value + intercept_;
    }

    double weighted_value = weightDatum(value, x_weight_);
    double eval = slope_ * weighted_value + intercept_;
    eval = unWeightDatum(eval, y_weight_);
    return eval;
  }

  void TransformationModelLinear::invert()
  {
    if (slope_ == 0)
    {
      throw Exception::DivisionByZero(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION);
    }
    intercept_ = -intercept_ / slope_;
    slope_ = 1.0 / slope_;

    // invert the weights:
    std::swap(x_datum_min_,y_datum_min_);
    std::swap(x_datum_max_,y_datum_max_);
    std::swap(x_weight_,y_weight_);

    // update parameters:
    params_.setValue("slope", slope_);
    params_.setValue("intercept", intercept_);
    params_.setValue("x_weight", x_weight_);
    params_.setValue("y_weight", y_weight_);
    params_.setValue("x_datum_min", x_datum_min_);
    params_.setValue("x_datum_max", x_datum_max_);
    params_.setValue("y_datum_min", y_datum_min_);
    params_.setValue("y_datum_max", y_datum_max_);
  }

  void TransformationModelLinear::getParameters(double& slope, double& intercept, String& x_weight, String& y_weight, double& x_datum_min, double& x_datum_max, double& y_datum_min, double& y_datum_max) const
  {
    slope = slope_;
    intercept = intercept_;
    x_weight = x_weight_;
    y_weight = y_weight_;
    x_datum_min = x_datum_min_;
    x_datum_max = x_datum_max_;
    y_datum_min = y_datum_min_;
    y_datum_max = y_datum_max_;
  }

  void TransformationModelLinear::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("symmetric_regression", "false", "Perform linear regression"
                                                     " on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.");
    params.setValidStrings("symmetric_regression",
                           {"true","false"});
    params.setValue("x_weight", "x", "Weight x values");
    params.setValidStrings("x_weight",
                           {"1/x","1/x2","ln(x)","x"});
    params.setValue("y_weight", "y", "Weight y values");
    params.setValidStrings("y_weight",
                           {"1/y","1/y2","ln(y)","y"});
    params.setValue("x_datum_min", 1e-15, "Minimum x value");
    params.setValue("x_datum_max", 1e15, "Maximum x value");
    params.setValue("y_datum_min", 1e-15, "Minimum y value");
    params.setValue("y_datum_max", 1e15, "Maximum y value");
  }

} // namespace
