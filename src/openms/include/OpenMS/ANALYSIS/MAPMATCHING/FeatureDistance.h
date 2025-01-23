// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Clemens Groepl, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/BaseFeature.h>

#include <limits>

namespace OpenMS
{
/**
   @brief A functor class for the calculation of distances between features or consensus features.

     It implements a customizable distance function of the following form:

     @f[
     w_{RT} \cdot \left( \frac{\left| RT_1 - RT_2 \right|}{\Delta RT_{max}} \right)^{p_{RT}} +
     w_{MZ} \cdot \left( \frac{\left| MZ_1 - MZ_2 \right|}{\Delta MZ_{max}} \right)^{p_{MZ}} +
     w_{int} \cdot \left( \frac{\left| int_1 - int_2 \right|}{int_{max}} \right)^{p_{int}}
     @f]

     This function returns a normalized distance between zero and one (unless constraints are violated, see below).

     @f$ RT_i @f$, @f$ MZ_i @f$, and @f$ int_i @f$ are the RT, m/z, and intensity values of the respective feature.

     Constraints are: @f$ {\Delta RT_{max}}, {\Delta MZ_{max}} @f$ and @f$ int_{max} @f$.
     If an absolute difference exceeds the specified maximum, the behavior depends on the value used for @p check_constraints in the constructor:
     If "false" (i.e., no constraints), the distance in that dimension may become greater than 1; if "true", @ref infinity is returned as overall distance.

     @f$ {\Delta RT_{max}} @f$ and @f$ {\Delta MZ_{max}} @f$ are the maximum allowed differences in RT and m/z, respectively.
     They are specified by the parameters @p distance_RT:max_difference and @p distance_MZ:max_difference, and are used for normalization,
     i.e., the observed RT or m/z differences of the feature pair are scaled relative to this value.

     @f$ int_{max} @f$ is the intensity which yields a normalized intensity of 1. This parameter is not settable via user params,
     but is set in the constructor (via parameter @p max_intensity), since it depends on the data at hand.

     @f$ p_X @f$ is the exponent for the distance in dimension X, specified by the parameter @p distance_X:exponent.
     Normalized differences (between (0, 1) unless unconstrained) are taken to this power. This makes it possible to compare values using linear, quadratic, etc. distance.

     @f$ w_X @f$ is the weight of final distance in dimension X, specified by the parameter @p distance_X:weight. The weights can be used to increase or decrease
     the contribution of RT, m/z, or intensity in the distance function.
     (The default weight for the intensity dimension is zero, i.e. intensity is not considered by default. However, @f$ int_{max} @f$ is still a constraint and
      should be set sensibly in the c'tor.)

     By default, two features are paired only if they have the same charge state (or at least one unknown charge '0') - otherwise, @ref infinity is returned.
     This behavior can be changed by the @p ignore_charge parameter.


     @note Peptide identifications annotated to features are not taken into account here,
           because they are stored in a format that is not suitable for rapid comparison.

   @htmlinclude OpenMS_FeatureDistance.parameters

   @ingroup FeatureGrouping
*/
  class OPENMS_DLLAPI FeatureDistance :
    public DefaultParamHandler
  {
public:
    /// Value to return if max. difference is exceeded or if charge states don't match
    static const double infinity;

    /**
       @brief Constructor

       @param max_intensity Maximum intensity of features (for normalization)
       @param force_constraints Check "max. difference" constraints given in the parameters and return @ref infinity if violated?
    */
    FeatureDistance(double max_intensity = 1.0,
                    bool force_constraints = false);

    /// Destructor
    ~FeatureDistance() override;

    /// Assignment operator
    FeatureDistance & operator=(const FeatureDistance & other);

    /**
       @brief Evaluation operator - checks constraints and computes the distance between two features

       @returns In the first element, whether constraints were satisfied; in
       the second element, the distance (@ref infinity if constraints were
       violated and @ref force_constraints_ is true).
    */
    std::pair<bool, double> operator()(const BaseFeature & left,
                                       const BaseFeature & right);

protected:

    /// Structure for storing distance parameters
    struct DistanceParams_
    {
      DistanceParams_() {}

      DistanceParams_(const String & what, const Param & global)
      {
        Param param = global.copy("distance_" + what + ":", true);
        if (what == "MZ")
        {
          max_diff_ppm = (param.getValue("unit") == "ppm");
        }
        else
        {
          max_diff_ppm = false;
        }

        max_difference = param.getValue("max_difference");
        exponent = param.getValue("exponent");
        weight = param.getValue("weight");
        norm_factor = 1 / max_difference;

        relevant = (weight != 0.0) && (exponent != 0.0);
        if (!relevant)
        {
            weight = 0.0;
        }
      }

      double max_difference, exponent, weight, norm_factor;
      bool max_diff_ppm, relevant;
    };

    /// Docu in base class
    void updateMembers_() override;

    /// Computes a distance component given absolute difference and parameters
    inline double distance_(double diff, const DistanceParams_ & params) const;

    /// Storage of parameters for the individual distance components
    DistanceParams_ params_rt_, params_mz_, params_intensity_;

    /// Reciprocal value of the total weight in the distance function
    double total_weight_reciprocal_;

    /// Maximum intensity of features (for normalization)
    double max_intensity_;

    /// Compute a distance even if charge states don't match?
    bool ignore_charge_;

    /// Compute a distance even if adducts don't match?
    bool ignore_adduct_;

    /// Always return @ref infinity if "max. difference" constraints are not met?
    bool force_constraints_;

    /// Log-transform intensities when computing intensity distance?
    bool log_transform_;
  };

} // namespace OpenMS

