// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser
// $Authors: Clemens Groepl, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace std;

namespace OpenMS
{

  const double FeatureDistance::infinity =
    std::numeric_limits<double>::infinity();


  FeatureDistance::FeatureDistance(double max_intensity,
                                   bool force_constraints) :
    DefaultParamHandler("FeatureDistance"),
    params_rt_(),
    params_mz_(),
    params_intensity_(),
    max_intensity_(max_intensity),
    force_constraints_(force_constraints),
    log_transform_(false)
  {
    defaults_.setValue("distance_RT:max_difference", 100.0, "Never pair features with a larger RT distance (in seconds).");
    defaults_.setMinFloat("distance_RT:max_difference", 0.0);
    defaults_.setValue("distance_RT:exponent", 1.0, "Normalized RT differences ([0-1], relative to 'max_difference') are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)", {"advanced"});
    defaults_.setMinFloat("distance_RT:exponent", 0.0);
    defaults_.setValue("distance_RT:weight", 1.0, "Final RT distances are weighted by this factor", {"advanced"});
    defaults_.setMinFloat("distance_RT:weight", 0.0);
    defaults_.setSectionDescription("distance_RT", "Distance component based on RT differences");

    defaults_.setValue("distance_MZ:max_difference", 0.3, "Never pair features with larger m/z distance (unit defined by 'unit')");
    defaults_.setMinFloat("distance_MZ:max_difference", 0.0);
    defaults_.setValue("distance_MZ:unit", "Da", "Unit of the 'max_difference' parameter");
    defaults_.setValidStrings("distance_MZ:unit", {"Da","ppm"});
    defaults_.setValue("distance_MZ:exponent", 2.0, "Normalized ([0-1], relative to 'max_difference') m/z differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)", {"advanced"});
    defaults_.setMinFloat("distance_MZ:exponent", 0.0);
    defaults_.setValue("distance_MZ:weight", 1.0, "Final m/z distances are weighted by this factor", {"advanced"});
    defaults_.setMinFloat("distance_MZ:weight", 0.0);
    defaults_.setSectionDescription("distance_MZ", "Distance component based on m/z differences");

    defaults_.setValue("distance_intensity:exponent", 1.0, "Differences in relative intensity ([0-1]) are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)", {"advanced"});
    defaults_.setMinFloat("distance_intensity:exponent", 0.0);
    defaults_.setValue("distance_intensity:weight", 0.0, "Final intensity distances are weighted by this factor", {"advanced"});
    defaults_.setMinFloat("distance_intensity:weight", 0.0);
    defaults_.setValue("distance_intensity:log_transform", "disabled", "Log-transform intensities? If disabled, d = |int_f2 - int_f1| / int_max. If enabled, d = |log(int_f2 + 1) - log(int_f1 + 1)| / log(int_max + 1))", {"advanced"});
    defaults_.setValidStrings("distance_intensity:log_transform", {"enabled","disabled"});
    defaults_.setSectionDescription("distance_intensity", "Distance component based on differences in relative intensity (usually relative to highest peak in the whole data set)");
    defaults_.setValue("ignore_charge", "false", "false [default]: pairing requires equal charge state (or at least one unknown charge '0'); true: Pairing irrespective of charge state");
    defaults_.setValidStrings("ignore_charge", {"true","false"});
    defaults_.setValue("ignore_adduct", "true", "true [default]: pairing requires equal adducts (or at least one without adduct annotation); true: Pairing irrespective of adducts");
    defaults_.setValidStrings("ignore_adduct", {"true","false"});


    defaultsToParam_();
  }

  FeatureDistance::~FeatureDistance() = default;

  FeatureDistance & FeatureDistance::operator=(const FeatureDistance & other)
  {
    DefaultParamHandler::operator=(other);

    max_intensity_ = other.max_intensity_;
    force_constraints_ = other.force_constraints_;
    updateMembers_();     // this sets all other member variables

    return *this;
  }

  void FeatureDistance::updateMembers_()
  {
    params_rt_ = DistanceParams_("RT", param_);
    params_mz_ = DistanceParams_("MZ", param_);

    log_transform_ = (param_.getValue("distance_intensity:log_transform") == "enabled");
    if (log_transform_)
    {
      // this parameter is not set by the user, but comes from the data:
      param_.setValue("distance_intensity:max_difference", Math::linear2log(max_intensity_));
    }
    else
    {
      // this parameter is not set by the user, but comes from the data:
      param_.setValue("distance_intensity:max_difference", max_intensity_);
    }
    params_intensity_ = DistanceParams_("intensity", param_);
    total_weight_reciprocal_ = 1 / (params_rt_.weight + params_mz_.weight +
                                    params_intensity_.weight);
    ignore_charge_ = param_.getValue("ignore_charge").toBool();
    ignore_adduct_ = param_.getValue("ignore_adduct").toBool();
  }

  double FeatureDistance::distance_(double diff, const DistanceParams_ & params) const
  {
    // manually querying for ^1 and ^2, since pow(x,2.0) is REALLY expensive and ^1 and ^2 are the defaults (so are likely to be used)
    if (params.exponent == 1)
    {
      return diff * params.norm_factor * params.weight;
    }
    else if (params.exponent == 2)
    {
      double tmp(diff * params.norm_factor);
      return tmp * tmp * params.weight;
    }
    else
    {
      // this pow() is REALLY expensive, since it uses a 'double' as exponent,
      // using 'int' will make it faster, but we will loose fractional
      // exponents (might be useful?).
      return pow(diff * params.norm_factor, params.exponent) * params.weight;
    }
  }

  pair<bool, double> FeatureDistance::operator()(const BaseFeature & left,
                                                 const BaseFeature & right)
  {
    if (!ignore_charge_)
    {
      Int charge_left = left.getCharge(), charge_right = right.getCharge();
      if (charge_left != charge_right)
      {
        if ((charge_left != 0) && (charge_right != 0))
        {
          return make_pair(false, infinity);
        }
      }
    }

    if (!ignore_adduct_)
    {
      if (left.metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS) && right.metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
      {
        if (EmpiricalFormula(left.getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)) != EmpiricalFormula(right.getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)))
        {
          return make_pair(false, infinity);
        }
      }
    }

    bool valid = true;

    // check m/z difference constraint:
    double left_mz = left.getMZ(), right_mz = right.getMZ();
    double dist_mz = fabs(left_mz - right_mz);
    double max_diff_mz = params_mz_.max_difference;
    if (params_mz_.max_diff_ppm) // compute absolute difference (in Da/Th)
    {
      max_diff_mz *= left_mz * 1e-6;
      // overwrite this parameter - it will be recomputed each time anyway:
      params_mz_.norm_factor = 1 / max_diff_mz;
    }

    if (dist_mz > max_diff_mz)
    {
      if (force_constraints_)
      {
        return make_pair(false, infinity);
      }
      valid = false;
    }

    // check RT difference constraint:
    double dist_rt = fabs(left.getRT() - right.getRT());
    if (dist_rt > params_rt_.max_difference)
    {
      if (force_constraints_)
      {
        return make_pair(false, infinity);
      }
      valid = false;
    }

    dist_rt = distance_(dist_rt, params_rt_);
    dist_mz = distance_(dist_mz, params_mz_);

    double dist_intensity = 0.0;
    if (params_intensity_.relevant)     // not by default, so worth checking
    {
      if (log_transform_)
      {
        dist_intensity = fabs(Math::linear2log(left.getIntensity()) - Math::linear2log(right.getIntensity()));
      }
      else
      {
        dist_intensity = fabs(left.getIntensity() - right.getIntensity());
      }

      dist_intensity = distance_(dist_intensity, params_intensity_);
    }

    double dist = dist_rt + dist_mz + dist_intensity;
    dist *= total_weight_reciprocal_;

    return make_pair(valid, dist);
  }

}
