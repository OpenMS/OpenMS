// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace std;

namespace OpenMS
{

  TransformationModelBSpline::TransformationModelBSpline(
    const TransformationModel::DataPoints& data, const Param& params) :
    spline_(nullptr)
  {
    // parameter handling/checking:
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    if (data.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "'b_spline' model requires more data");
    }
    Size boundary_condition = params_.getValue("boundary_condition");

    BSpline2d::BoundaryCondition bound_cond = 
      static_cast<BSpline2d::BoundaryCondition>(boundary_condition);
    vector<double> x(data.size()), y(data.size());
    xmin_ = data[0].first;
    xmax_ = xmin_;
    for (Size i = 0; i < data.size(); ++i)
    {
      x[i] = data[i].first;
      y[i] = data[i].second;
      if (x[i] < xmin_) xmin_ = x[i];
      else if (x[i] > xmax_) xmax_ = x[i];
    }
    double wavelength = params_.getValue("wavelength");
    if (wavelength > (xmax_ - xmin_))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "B-spline 'wavelength' can't be larger than the data range (here: " + String(xmax_ - xmin_) + ").", String(wavelength));
    }

    // since we can't initialize a BSpline2d object in the init list (no c'tor
    // that doesn't require preparation of data), we have to use a pointer:
    spline_ = new BSpline2d(x, y, wavelength, bound_cond, 
                            params_.getValue("num_nodes"));

    if (!spline_->ok())
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                   "TransformationModelBSpline", 
                                   "Unable to fit B-spline to data points.");
    }

    // extrapolation:
    std::string extrapolate = params_.getValue("extrapolate");
    if (extrapolate == "b_spline")
    {
      extrapolate_ = EX_BSPLINE;
    }
    else if (extrapolate == "global_linear")
    {
      extrapolate_ = EX_GLOBAL_LINEAR;
      TransformationModelLinear lm(data, Param());
      String x_weight, y_weight;
      double x_datum_min, x_datum_max, y_datum_min, y_datum_max;
      lm.getParameters(slope_min_, offset_min_, x_weight, y_weight, x_datum_min, x_datum_max, y_datum_min, y_datum_max);
      slope_max_ = slope_min_;
      // extrapolation (left/right) considers xmin_/xmax_ as the origin (x = 0):
      offset_min_ = lm.evaluate(xmin_);
      offset_max_ = lm.evaluate(xmax_);
    }
    else // "linear" or "constant"
    {
      offset_min_ = spline_->eval(xmin_);
      offset_max_ = spline_->eval(xmax_);
      if (extrapolate == "constant") 
      {
        extrapolate_ = EX_CONSTANT;
      }
      else // "linear"
      {
        extrapolate_ = EX_LINEAR;
        slope_min_ = spline_->derivative(xmin_);
        slope_max_ = spline_->derivative(xmax_);
      }
    }
  }


  TransformationModelBSpline::~TransformationModelBSpline()
  {
    if (spline_) delete spline_;
  }

  double TransformationModelBSpline::evaluate(double value) const
  {
    if ((value < xmin_) && (extrapolate_ != EX_BSPLINE)) // extrapolate (left)
    {
      if (extrapolate_ == EX_CONSTANT)
      {
        return offset_min_;
      }
      else // "EX_LINEAR" or "EX_GLOBAL_LINEAR"
      {
        return offset_min_ - slope_min_ * (xmin_ - value);
      }
    }
    if ((value > xmax_) && (extrapolate_ != EX_BSPLINE)) // extrapolate (right)
    {
      if (extrapolate_ == EX_CONSTANT)
      {
        return offset_max_;
      }
      else // "EX_LINEAR" or "EX_GLOBAL_LINEAR"
      {
        return offset_max_ + slope_max_ * (value - xmax_);
      }
    }
    return spline_->eval(value);
  }

  void TransformationModelBSpline::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("wavelength", 0.0, "Determines the amount of smoothing by setting the number of nodes for the B-spline. The number is chosen so that the spline approximates a low-pass filter with this cutoff wavelength. The wavelength is given in the same units as the data; a higher value means more smoothing. '0' sets the number of nodes to twice the number of input points.");
    params.setMinFloat("wavelength", 0.0);
    params.setValue("num_nodes", 5, "Number of nodes for B-spline fitting. Overrides 'wavelength' if set (to two or greater). A lower value means more smoothing.");
    params.setMinInt("num_nodes", 0);
    params.setValue("extrapolate", "linear", "Method to use for extrapolation beyond the original data range. 'linear': Linear extrapolation using the slope of the B-spline at the corresponding endpoint. 'b_spline': Use the B-spline (as for interpolation). 'constant': Use the constant value of the B-spline at the corresponding endpoint. 'global_linear': Use a linear fit through the data (which will most probably introduce discontinuities at the ends of the data range).");
    params.setValidStrings("extrapolate", {"linear","b_spline","constant","global_linear"});
    params.setValue("boundary_condition", 2, "Boundary condition at B-spline endpoints: 0 (value zero), 1 (first derivative zero) or 2 (second derivative zero)");
    params.setMinInt("boundary_condition", 0);
    params.setMaxInt("boundary_condition", 2);
  }

} // namespace
