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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>

namespace OpenMS
{
	/**
	@brief Generic description of a coordinate transformation.
		
	This description primarily stores data points (coordinate pairs) from which a @ref TransformationModel "transformation model" can be estimated. Applying the transformation to a coordinate (via @p apply) then means evaluating the model at that coordinate.

	The following models are available:
	- @p none (TransformationModel): \f$ f(x) = x \f$ (identity)
	- @p identity: Same as @p none, but intended for reference files (used to indicate that no other model should be fit, because the identity is already optimal).
	- @p linear (TransformationModelLinear): \f$ f(x) = slope * x + intercept \f$
	- @p interpolated (TransformationModelInterpolated): Interpolation between pairs, extrapolation using first and last pair. Supports different interpolation types.
	- @p b_spline (TransformationModelBSpline): Smoothing cubic B-spline.

	@remark TransformationDescription stores data points, TransformationModel stores parameters. That way, data can be modeled using different models/parameters, and models can still keep a representation of the data in the format they need (if at all).
		
	@ingroup MapAlignment
	*/
	class OPENMS_DLLAPI TransformationDescription
	{
		// friend class MapAlignmentAlgorithm;

	 public:
			
		/// Coordinate pair
		typedef TransformationModel::DataPoint DataPoint;
		/// Vector of coordinate pairs
		typedef TransformationModel::DataPoints DataPoints;
			
		/// Default constructor
		TransformationDescription();
		/// Constructor from data
		TransformationDescription(const DataPoints& data);
		/// Destructor
		~TransformationDescription();
				
		/// Copy constructor 
		TransformationDescription(const TransformationDescription& rhs);
		/// Assignment operator
		TransformationDescription& operator=(const TransformationDescription& rhs);
				
		/// Fits a model to the data
		void fitModel(const String& model_type, const Param& params=Param());

		/**
			 @brief Applies the transformation to @p value.

			 Returns the result of evaluating the fitted model at @p value.
			 Returns @p value unchanged if no model was fitted.
		*/
		DoubleReal apply(DoubleReal value) const;

		/// Gets the type of the fitted model
		const String& getModelType() const;

		/// Gets the possible types of models
		static void getModelTypes(StringList& result);

		/**
			 @brief Sets the data points
			 
			 Removes the model that was previously fitted to the data (if any).
		*/
		void setDataPoints(const DataPoints& data);

		/// Returns the data points
		const DataPoints& getDataPoints() const;

		/// Non-mutable access to the model parameters
		void getModelParameters(Param& params) const;
			
		/// Computes an (approximate) inverse of the transformation
		void invert();

	protected:
		/// Data points
		DataPoints data_;
		/// Type of model
		String model_type_;
		/// Pointer to model
		TransformationModel* model_;
	};
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONDESCRIPTION_H

