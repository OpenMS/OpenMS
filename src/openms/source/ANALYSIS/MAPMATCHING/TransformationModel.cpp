// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::find
#include <cmath>    // std::log

namespace OpenMS
{

  TransformationModel::TransformationModel(const TransformationModel::DataPoints&, const Param& p) :
    params_(p),
    x_weight_("x"),
    x_datum_min_(0),
    x_datum_max_(0),
    y_weight_("y"),
    y_datum_min_(0),
    y_datum_max_(0),
    weighting_(false)
  {
    // get x datum ranges
    x_datum_min_ = params_.exists("x_datum_min") ? (double)params_.getValue("x_datum_min") : 1e-15;
    x_datum_max_ = params_.exists("x_datum_max") ? (double)params_.getValue("x_datum_max") : 1e15;

    // get y datum ranges
    y_datum_min_ = params_.exists("y_datum_min") ? (double)params_.getValue("y_datum_min") : 1e-15;
    y_datum_max_ = params_.exists("y_datum_max") ? (double)params_.getValue("y_datum_max") : 1e15;

    // TrafoXML's prior to OpenMS 3.0 have x/y_weight = "" if unweighted
    x_weight_ = params_.exists("x_weight") && (params_.getValue("x_weight") != "") ? String(params_.getValue("x_weight").toString()) : "x";
    y_weight_ = params_.exists("y_weight") && (params_.getValue("y_weight") != "") ? String(params_.getValue("y_weight").toString()) : "y";

    std::vector<String> valid_x_weights = getValidXWeights();
    std::vector<String> valid_y_weights = getValidYWeights();
    if (x_weight_ != "x" && !checkValidWeight(x_weight_, valid_x_weights))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value '" + x_weight_ + "' is not a valid weight parameter for x values.");
    }
    if (y_weight_ != "y" && !checkValidWeight(y_weight_, valid_y_weights))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value '" + y_weight_ + "' is not a valid weight parameter for y values.");
    }

    // easily remember whether we do weighting or not
    weighting_ = !(x_weight_ == "x" && y_weight_ == "y"); 
  }

  TransformationModel::~TransformationModel() = default;

  double TransformationModel::evaluate(double value) const
  {
    return value;
  }

  const Param& TransformationModel::getParameters() const
  {
    return params_;
  }

  void TransformationModel::getDefaultParameters(Param& params)
  {
    params.clear();
  }
  
  void TransformationModel::weightData(TransformationModel::DataPoints& data)
  {
    if (!weighting_ ) return;

    // weight x values 
    if (x_weight_ != "x" && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        // check x datum ranges
        data[i].first = checkDatumRange(data[i].first,x_datum_min_,x_datum_max_);
        // weight x datum
        data[i].first = weightDatum(data[i].first, x_weight_);
      }
    }

    // weight y values
    if (y_weight_ != "y" && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        // check y datum ranges
        data[i].second = checkDatumRange(data[i].second,y_datum_min_, y_datum_max_);
        // weight y datum
        data[i].second = weightDatum(data[i].second, y_weight_);
      }
    } 
  }
  
  void TransformationModel::unWeightData(TransformationModel::DataPoints& data)
  {
    if (!weighting_ ) return;

    // unweight x values 
    if (x_weight_ != "x" && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        data[i].first = unWeightDatum(data[i].first, x_weight_);
      }
    }
    // unweight y values
    if (y_weight_ != "y" && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        data[i].second = unWeightDatum(data[i].second, y_weight_);
      }
    }  
  }

  bool TransformationModel::checkValidWeight(const String& weight, const std::vector<String>& valid_weights) const
  {    
    bool valid = false;
    if (std::find(valid_weights.begin(), valid_weights.end(), weight) != valid_weights.end())
    {
      valid=true;
    }
    else
    {
      OPENMS_LOG_INFO << "weight " + weight + " is not supported.";
    }
    return valid;
  }

  double TransformationModel::checkDatumRange(const double& datum, const double& datum_min, const double& datum_max)
  {    
    double datum_checked = datum;
    if (datum >= datum_max)
    {
      OPENMS_LOG_INFO << "datum " << datum << " is out of range.";
      OPENMS_LOG_INFO << "datum will be truncated to " << datum_max << ".";
      datum_checked = datum_max;
    }
    else if (datum <= datum_min)
    {
      OPENMS_LOG_INFO << "datum " << datum << " is out of range.";
      OPENMS_LOG_INFO << "datum will be truncated to " << datum_min << ".";
      datum_checked = datum_min;
    }
    return datum_checked;
  }
  
  std::vector<String> TransformationModel::getValidXWeights() const
  {
    std::vector<String> valid_weights{"1/x","1/x2","ln(x)","x"}; // == 1 disables weights
    return valid_weights;
  }
  
  std::vector<String> TransformationModel::getValidYWeights() const
  {
    std::vector<String> valid_weights{"1/y","1/y2","ln(y)","y"}; // == 1 disables weights
    return valid_weights;
  }

  double TransformationModel::weightDatum(const double& datum, const String& weight) const
  { 
    double datum_weighted = 0;   
    if (weight == "ln(x)")
    {
      datum_weighted = std::log(datum);
    }
    else if (weight == "ln(y)")
    {
      datum_weighted = std::log(datum);
    }
    else if (weight == "1/x")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/y")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/x2")
    {
      datum_weighted = 1/std::pow(datum,2);
    }
    else if (weight == "1/y2")
    {
      datum_weighted = 1/std::pow(datum,2);
    }
    else if (weight == "x" || weight == "y" )
    {
      datum_weighted = datum;
    }
    else
    {
      datum_weighted = datum;
      OPENMS_LOG_INFO << "weight " + weight + " not supported.";
      OPENMS_LOG_INFO << "no weighting will be applied.";
    }
    return datum_weighted;
  } 

  double TransformationModel::unWeightDatum(const double& datum, const String& weight) const
  { 
    double datum_weighted = 0;   
    if (weight == "ln(x)")
    {
      datum_weighted = std::exp(datum);
    }
    else if (weight == "ln(y)")
    {
      datum_weighted = std::exp(datum);
    }
    else if (weight == "1/x")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/y")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/x2")
    {
      datum_weighted = std::sqrt(1/std::abs(datum));
    }
    else if (weight == "1/y2")
    {
      datum_weighted = std::sqrt(1/std::abs(datum));
    }
    else if (weight == "x" || weight == "y")
    {
      datum_weighted = datum;
    }
    else
    {
      datum_weighted = datum;
      OPENMS_LOG_INFO << "weight " + weight + " not supported.";
      OPENMS_LOG_INFO << "no weighting will be applied.";
    }
    return datum_weighted;
  }

}
