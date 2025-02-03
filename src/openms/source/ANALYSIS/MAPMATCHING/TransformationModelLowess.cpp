// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/PROCESSING/SMOOTHING/FastLowessSmoothing.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  bool cmpFirstDimension(const TransformationModel::DataPoint& x, const TransformationModel::DataPoint& y)
  {
    return (x.first < y.first);
  }

  TransformationModelLowess::TransformationModelLowess(
      const TransformationModel::DataPoints& data_,
      const Param& params) : model_(nullptr)
  {
    // parameter handling/checking:
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    if (data_.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "'lowess' model requires more data");
    }

    // TODO copy ... 
    TransformationModel::DataPoints data(data_);

    // sort data
    std::sort(data.begin(), data.end(), cmpFirstDimension);

    vector<double> x(data.size()), y(data.size()), result(data.size());
    double xmin_ = data[0].first;
    double xmax_ = xmin_;
    for (Size i = 0; i < data.size(); ++i)
    {
      x[i] = data[i].first;
      y[i] = data[i].second;
      if (x[i] < xmin_) 
      {
        xmin_ = x[i];
      }
      else if (x[i] > xmax_)
      {
        xmax_ = x[i];
      }
    }

    double span = params_.getValue("span");
    int nsteps = params_.getValue("num_iterations");
    double delta = params_.getValue("delta");
    
    if (delta < 0.0)
    {
      delta = (xmax_ - xmin_) * 0.01; // automatically determine delta
    }

    FastLowessSmoothing::lowess(x, y, span, nsteps, delta, result);

    TransformationModel::DataPoints data_out;
    for (Size i = 0; i < result.size(); ++i)
    {
      data_out.push_back( std::make_pair(x[i], result[i]) );
    }

    // TODO thin out data here ? we may not need that many points here to interpolate ...  it is enough if we store a few datapoints

    Param p;
    TransformationModelInterpolated::getDefaultParameters(p);
    /// p.setValue("interpolation_type", "cspline"); // linear interpolation between lowess pts
    /// p.setValue("extrapolation_type", "four-point-linear");
    p.setValue("interpolation_type", params_.getValue("interpolation_type"));
    p.setValue("extrapolation_type", params_.getValue("extrapolation_type"));

    // create new interpolation model based on the lowess data
    model_ = new TransformationModelInterpolated(data_out, p);
  }

  TransformationModelLowess::~TransformationModelLowess()
  {
    if (model_) delete model_;
  }

  void TransformationModelLowess::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("span", 2/3.0, "Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.");
    params.setMinFloat("span", 0.0);
    params.setMaxFloat("span", 1.0);

    params.setValue("num_iterations", 3, "Number of robustifying iterations for lowess fitting.");
    params.setMinInt("num_iterations", 0);

    params.setValue("delta", -1.0, "Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.");

    params.setValue("interpolation_type", "cspline", "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation");
    params.setValidStrings("interpolation_type", {"linear","cspline","akima"});

    params.setValue("extrapolation_type", "four-point-linear", "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation.");
    params.setValidStrings("extrapolation_type", {"two-point-linear","four-point-linear","global-linear"});
  }

}
