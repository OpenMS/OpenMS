// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREDISTANCE_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREDISTANCE_H

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

	 @f$ {\Delta RT_{max}} @f$ and @f$ {\Delta MZ_{max}} @f$ are the maximum allowed differences in RT and m/z, respectively. They are specified by the parameters @p distance_RT:max_difference and @p distance_MZ:max_difference, and are used for normalization. If an absolute difference exceeds the specified maximum, the behavior depends on the value used for @p check_constraints in the constructor: If "false", the distance will be computed normally, but may become greater than 1; if "true", the fixed value @ref infinity is returned.

	 @f$ int_{max} @f$ is the maximum intensity that can occur for features compared by this distance function. It is not a parameter specified by the user, but depends on the data at hand and is thus set in the constructor (via parameter @p max_intensity).

	 @f$ w_X @f$ is the weight of distance component X, specified by the parameter @p distance_X:weight. The weights can be used to increase or decrease the contribution of RT, m/z, or intensity in the distance function. (Note that the default weight for the intensity component is zero, i.e. intensity is not considered by default.)

	 @f$ p_X @f$ is the exponent for distance component X, specified by the parameter @p distance_X:exponent. Normalized differences are taken to this power. This makes it possible to compare values using linear, quadratic, etc. distance.

	 By default, two features are only compared if they have the same charge state (or charge state 0 for "undefined") - otherwise, @ref infinity is returned. This behavior can be changed by the @p ignore_charge parameter.

	 @note Peptide identifications annotated to features are not taken into account here, because they are stored in a format that is not suitable for rapid comparison.

   @htmlinclude OpenMS_FeatureDistance.parameters

   @ingroup FeatureGrouping
*/
	class OPENMS_DLLAPI FeatureDistance: public DefaultParamHandler
	{
	public:
		/// Value to return if max. difference is exceeded or if charge states don't match
		static const DoubleReal infinity;

		/**
			 @brief Constructor
			 
			 @param max_intensity Maximum intensity of features (for normalization)
			 @param check_constraints Check "max. difference" constraints given in the parameters and return @ref infinity if violated?
		*/
		FeatureDistance(DoubleReal max_intensity = 1.0,
										bool force_constraints = false);
		
		/// Desctructor
		virtual ~FeatureDistance();

		/// Assignment operator
		FeatureDistance& operator=(const FeatureDistance& other);

		/**
			 @brief Evaluation operator - checks constraints and computes the distance between two features

			 @returns In the first element, whether constraints were satisfied; in the second element, the distance (@ref infinity if constraints were violated and @ref force_constraints_ is true).
		*/
		std::pair<bool, DoubleReal> operator()(const BaseFeature& left, 
																					 const BaseFeature& right);

	protected:
		/// Structure for storing distance parameters
		struct DistanceParams_;

		/// Docu in base class
		void updateMembers_();

		/// Computes a distance component given absolute difference and parameters
		inline DoubleReal distance_(DoubleReal diff, const DistanceParams_* params);

		/// Storage of parameters for the individual distance components
		DistanceParams_ *params_rt_, *params_mz_, *params_intensity_;
		
		/// Reciprocal value of the total weight in the distance function
		DoubleReal total_weight_reciprocal_;

		/// Maximum intensity of features (for normalization)
		DoubleReal max_intensity_;

		/// Compute a distance even if charge states don't match?
		bool ignore_charge_;

		/// Always return @ref infinity if "max. difference" constraints are not met?
		bool force_constraints_;
	};

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_FEATUREDISTANCE_H
