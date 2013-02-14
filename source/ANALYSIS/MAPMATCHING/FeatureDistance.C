// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Hendrik Weisser
// $Authors: Clemens Groepl, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

using namespace std;

namespace OpenMS
{

  const DoubleReal FeatureDistance::infinity =
    std::numeric_limits<DoubleReal>::infinity();


  FeatureDistance::FeatureDistance(DoubleReal max_intensity,
                                   bool force_constraints) :
    DefaultParamHandler("FeatureDistance"),
    params_rt_(), params_mz_(), params_intensity_(),
    max_intensity_(max_intensity), force_constraints_(force_constraints)
  {
    defaults_.setValue("distance_RT:max_difference", 100.0, "Maximum allowed difference in RT in seconds");
    defaults_.setMinFloat("distance_RT:max_difference", 0.0);
    defaults_.setValue("distance_RT:exponent", 1.0, "Normalized RT differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)", StringList::create("advanced"));
    defaults_.setMinFloat("distance_RT:exponent", 0.0);
    defaults_.setValue("distance_RT:weight", 1.0, "RT distances are weighted by this factor", StringList::create("advanced"));
    defaults_.setMinFloat("distance_RT:weight", 0.0);
    defaults_.setSectionDescription("distance_RT", "Distance component based on RT differences");

    defaults_.setValue("distance_MZ:max_difference", 0.3, "Maximum allowed difference in m/z (unit defined by 'unit')");
    defaults_.setMinFloat("distance_MZ:max_difference", 0.0);
    defaults_.setValue("distance_MZ:unit", "Da", "Unit of the 'max_difference' parameter");
    defaults_.setValidStrings("distance_MZ:unit", StringList::create("Da,ppm"));
    defaults_.setValue("distance_MZ:exponent", 2.0, "Normalized m/z differences are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)", StringList::create("advanced"));
    defaults_.setMinFloat("distance_MZ:exponent", 0.0);
    defaults_.setValue("distance_MZ:weight", 1.0, "m/z distances are weighted by this factor", StringList::create("advanced"));
    defaults_.setMinFloat("distance_MZ:weight", 0.0);
    defaults_.setSectionDescription("distance_MZ", "Distance component based on m/z differences");

    defaults_.setValue("distance_intensity:exponent", 1.0, "Differences in relative intensity are raised to this power (using 1 or 2 will be fast, everything else is REALLY slow)", StringList::create("advanced"));
    defaults_.setMinFloat("distance_intensity:exponent", 0.0);
    defaults_.setValue("distance_intensity:weight", 0.0, "Distances based on relative intensity are weighted by this factor", StringList::create("advanced"));
    defaults_.setMinFloat("distance_intensity:weight", 0.0);
    defaults_.setSectionDescription("distance_intensity", "Distance component based on differences in relative intensity");

    defaults_.setValue("ignore_charge", "false", "Compare features normally even if their charge states are different");
    defaults_.setValidStrings("ignore_charge", StringList::create("true,false"));

    defaultsToParam_();
  }

  FeatureDistance::~FeatureDistance()
  {
  }

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
    // this parameter is not set by the user, but comes from the data:
    param_.setValue("distance_intensity:max_difference", max_intensity_);
    params_intensity_ = DistanceParams_("intensity", param_);
    total_weight_reciprocal_ = 1 / (params_rt_.weight + params_mz_.weight +
                                    params_intensity_.weight);
    ignore_charge_ = String(param_.getValue("ignore_charge")) == "true";
  }

  DoubleReal FeatureDistance::distance_(DoubleReal diff, const DistanceParams_ & params) const
  {
    // manually querying for ^1 and ^2, since pow(x,2.0) is REALLY expensive and ^1 and ^2 are the defaults (so are likely to be used)
    if (params.exponent == 1)
      return diff * params.norm_factor * params.weight;
    else if (params.exponent == 2)
    {
      DoubleReal tmp(diff * params.norm_factor);
      return tmp * tmp * params.weight;
    }
    else // this pow() is REALLY expensive, since it uses a 'double' as exponent, using 'int' will make it faster,
    { // but we will loose fractional exponents (might be useful?)
      return pow(diff * params.norm_factor, params.exponent) * params.weight;
    }
  }

  pair<bool, DoubleReal> FeatureDistance::operator()(const BaseFeature & left,
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

    bool valid = true;

    // check m/z difference constraint:
    DoubleReal left_mz = left.getMZ(), right_mz = right.getMZ();
    DoubleReal dist_mz = fabs(left_mz - right_mz);
    DoubleReal max_diff_mz = params_mz_.max_difference;
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
    DoubleReal dist_rt = fabs(left.getRT() - right.getRT());
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

    DoubleReal dist_intensity = 0.0;
    if (params_intensity_.relevant)     // not by default, so worth checking
    {
      dist_intensity = fabs(left.getIntensity() - right.getIntensity());
      dist_intensity = distance_(dist_intensity, params_intensity_);
    }

    DoubleReal dist = dist_rt + dist_mz + dist_intensity;
    dist *= total_weight_reciprocal_;

    return make_pair(valid, dist);
  }

}
