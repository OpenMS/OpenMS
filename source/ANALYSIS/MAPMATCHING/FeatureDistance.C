// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

using namespace std;

namespace OpenMS
{

	const DoubleReal FeatureDistance::infinity = 
		std::numeric_limits<DoubleReal>::infinity();


	struct FeatureDistance::DistanceParams_
	{
		DistanceParams_(const String& what, const Param& global)
		{
			Param param = global.copy("distance_" + what + ":", true);
			if (what == "MZ") max_diff_ppm = param.getValue("unit") == "ppm";
			else max_diff_ppm = false;
			max_difference = param.getValue("max_difference");
			exponent = param.getValue("exponent");
			weight = param.getValue("weight");
			norm_factor = 1 / max_difference;
			relevant = (weight != 0.0) && (exponent != 0.0);
			if (!relevant) weight = 0.0;
		}

		DoubleReal max_difference, exponent, weight, norm_factor;
		bool max_diff_ppm, relevant;
	};


	FeatureDistance::FeatureDistance(DoubleReal max_intensity, 
																	 bool force_constraints):
		DefaultParamHandler("FeatureDistance"), 
		params_rt_(0), params_mz_(0), params_intensity_(0),
		max_intensity_(max_intensity), force_constraints_(force_constraints)
	{
    defaults_.setValue("distance_RT:max_difference", 100.0, "Maximum allowed difference in RT");
    defaults_.setMinFloat("distance_RT:max_difference", 0.0);
		defaults_.setValue("distance_RT:exponent", 1.0, "Normalized RT differences are raised to this power", StringList::create("advanced"));
		defaults_.setMinFloat("distance_RT:exponent", 0.0);
		defaults_.setValue("distance_RT:weight", 1.0, "RT distances are weighted by this factor", StringList::create("advanced"));
		defaults_.setMinFloat("distance_RT:weight", 0.0);
		defaults_.setSectionDescription("distance_RT", "Distance component based on RT differences");

    defaults_.setValue("distance_MZ:max_difference", 0.3, "Maximum allowed difference in m/z (unit defined by 'mz_unit')");
    defaults_.setMinFloat("distance_MZ:max_difference", 0.0);
    defaults_.setValue("distance_MZ:unit", "Da", "Unit of the 'max_difference' parameter");
    defaults_.setValidStrings("distance_MZ:unit", StringList::create("Da,ppm"));
		defaults_.setValue("distance_MZ:exponent", 2.0, "Normalized m/z differences are raised to this power", StringList::create("advanced"));
		defaults_.setMinFloat("distance_MZ:exponent", 0.0);
		defaults_.setValue("distance_MZ:weight", 1.0, "m/z distances are weighted by this factor", StringList::create("advanced"));
		defaults_.setMinFloat("distance_MZ:weight", 0.0);
		defaults_.setSectionDescription("distance_MZ", "Distance component based on m/z differences");

		defaults_.setValue("distance_intensity:exponent", 1.0, "Differences in relative intensity are raised to this power", StringList::create("advanced"));
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
		delete params_rt_;
		delete params_mz_;
		delete params_intensity_;
	}


	FeatureDistance& FeatureDistance::operator=(const FeatureDistance& other)
	{
		DefaultParamHandler::operator=(other);

		max_intensity_ = other.max_intensity_;
		force_constraints_ = other.force_constraints_;
		updateMembers_(); // this sets all other member variables

		return *this;
	}


	void FeatureDistance::updateMembers_()
	{
		delete params_rt_;
		delete params_mz_;
		delete params_intensity_;
		params_rt_ = new DistanceParams_("RT", param_);
		params_mz_ = new DistanceParams_("MZ", param_);
		// this parameter is not set by the user, but comes from the data:
		param_.setValue("distance_intensity:max_difference", max_intensity_);
		params_intensity_ = new DistanceParams_("intensity", param_);
		total_weight_reciprocal_ = 1 / (params_rt_->weight + params_mz_->weight + 
																		params_intensity_->weight);
		ignore_charge_ = String(param_.getValue("ignore_charge")) == "true";
	}


	DoubleReal FeatureDistance::distance_(DoubleReal diff, 
																				const DistanceParams_* params)
	{
		return pow(diff * params->norm_factor, params->exponent) * params->weight;
	}


	pair<bool, DoubleReal> FeatureDistance::operator()(const BaseFeature& left, 
																										 const BaseFeature& right)
	{
		if (!ignore_charge_)
		{
			Int charge_left = left.getCharge(), charge_right = right.getCharge();
			if ((charge_left != 0) && (charge_right != 0) && 
					(charge_left != charge_right))
			{
				return make_pair(false, infinity);
			}
		}

		bool valid = true;
		DoubleReal dist_rt = 0.0, dist_mz = 0.0, dist_intensity = 0.0;
		
		// check RT difference constraint:
		dist_rt = abs(left.getRT() - right.getRT());
		if (dist_rt > params_rt_->max_difference)
		{
			if (force_constraints_) return make_pair(false, infinity);
			valid = false;
		}

		// check m/z difference constraint:
		DoubleReal left_mz = left.getMZ(), right_mz = right.getMZ();
		dist_mz = abs(left_mz - right_mz);
		DoubleReal max_diff_mz = params_mz_->max_difference;
		if (params_mz_->max_diff_ppm) // compute absolute difference (in Da/Th)
		{
			max_diff_mz *= left_mz * 1e-6;
			// overwrite this parameter - it will be recomputed each time anyway:
			params_mz_->norm_factor = 1 / max_diff_mz;
		}
		if (dist_mz > max_diff_mz)
		{
			if (force_constraints_) return make_pair(false, infinity);
			valid = false;
		}
		
		dist_rt = distance_(dist_rt, params_rt_);

		dist_mz = distance_(dist_mz, params_mz_);

		if (params_intensity_->relevant) // not by default, so worth checking
		{
			dist_intensity = abs(left.getIntensity() - right.getIntensity());
			dist_intensity = distance_(dist_intensity, params_intensity_);
		}
		
		DoubleReal dist = dist_rt + dist_mz + dist_intensity;
		dist *= total_weight_reciprocal_;

		return make_pair(valid, dist);
	}

}
