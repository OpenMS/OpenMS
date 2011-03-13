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
// $Maintainer: $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H
#define OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_interp.h>

namespace OpenMS
{
	/**
		 @brief Base class for transformation models

		 Implements the identity (no transformation). Parameters and data are ignored.

		 @ingroup MapAlignment
	*/
	class OPENMS_DLLAPI TransformationModel
	{
	public:
		/// Coordinate pair
		typedef std::pair<DoubleReal, DoubleReal> DataPoint;
		/// Vector of coordinate pairs
		typedef std::vector<DataPoint> DataPoints;
			
		/// Constructor
		TransformationModel() {};

		/// Alternative constructor (derived classes should implement this one!)
		TransformationModel(const TransformationModel::DataPoints&, 
												const Param&) {};
		
		/// Destructor
		virtual ~TransformationModel() {};

		/// Evaluates the model at the given value
		virtual DoubleReal evaluate(const DoubleReal value) const
		{
			return value;
		};

		/// Gets the parameters
		virtual void getParameters(Param& params) const
		{
			params.clear();
		};
	};


	/**
		 @brief Linear model for transformations

		 The model can be inferred from data or specified using explicit parameters. If data is given, a least squares fit is used to find the model parameters (slope and intercept). Depending on parameter @p symmetric_regression, a normal regression (@e y on @e x) or symmetric regression (@f$ y - x @f$ on @f$ y + x @f$) is performed.

		 Without data, the model can be specified by giving the parameters @p slope and @p intercept explicitly.

		 @ingroup MapAlignment
	*/
	class OPENMS_DLLAPI TransformationModelLinear: public TransformationModel
	{
	public:
		/**
			@brief Constructor
			
			@exception IllegalArgument is thrown if neither data points nor explicit parameters (slope/intercept) are given.
		*/
		TransformationModelLinear(const DataPoints& data, const Param& params);

		/// Destructor
		~TransformationModelLinear();

		/// Evaluates the model at the given value
		DoubleReal evaluate(const DoubleReal value) const;

		/// Gets the parameters (only if set explicitly, not if estimated)
		void getParameters(Param& params) const;

		/// Gets the "real" parameters
		void getParameters(DoubleReal& slope, DoubleReal& intercept) const;

		/**
			 @brief Computes the inverse

			 @exception DivisionByZero is thrown if the slope is zero.
		*/
		void invert();
		
	protected:
		/// Parameters of the linear model
		DoubleReal slope_, intercept_;
		/// Was the model estimated from data?
		bool data_given_;
		/// Use symmetric regression?
		bool symmetric_;
	};


	/**
		 @brief Interpolation model for transformations

		 Between the data points, the interpolation uses the neighboring points. Outside the range spanned by the points, we extrapolate using a line through the first and the last point.

		 Different types of interpolation (controlled by the parameter @p interpolation_type) are supported: "linear", "polynomial", "cspline", and "akima". Note that the number of required data points may differ between types.

		 @ingroup MapAlignment
	*/
	class OPENMS_DLLAPI TransformationModelInterpolated: 
		public TransformationModel
	{
	public:
		/**
			 @brief Constructor

			 @exception IllegalArgument is thrown if there are not enough data points or if an unknown interpolation type is given.
		*/
		TransformationModelInterpolated(const DataPoints& data, 
																		const Param& params);
		
		/// Destructor
		~TransformationModelInterpolated();

		/// Evaluates the model at the given value
		DoubleReal evaluate(const DoubleReal value) const;

		/// Gets the parameters
		void getParameters(Param& params) const;

		/// Gets allowed values for the parameter "interpolation_type"
		static void getInterpolationTypes(StringList& result);

	protected:
		/// Data coordinates
		std::vector<double> x_, y_;
		/// Number of data points
		size_t size_;
		/// Look-up accelerator
		gsl_interp_accel *acc_;
		/// Interpolation function
		gsl_interp *interp_;
		/// Linear model for extrapolation
		TransformationModelLinear* lm_;
	};


	/**
		 @brief B-spline model for transformations

		 In the range of the data points, the transformation is evaluated from a cubic smoothing spline fit to the points. The number of breakpoints is given as a parameter (@p num_breakpoints). Outside of this range, linear extrapolation through the last point with the slope of the spline at that point is used.

		 @ingroup MapAlignment
	*/
	class OPENMS_DLLAPI TransformationModelBSpline: public TransformationModel
	{
	public:
		/**
			 @brief Constructor

			 @exception IllegalArgument is thrown if not enough data points are given or if the required parameter @p num_breakpoints is missing.
		*/
		TransformationModelBSpline(const DataPoints& data, const Param& params);
		
		/// Destructor
		~TransformationModelBSpline();

		/// Evaluates the model at the given value
		DoubleReal evaluate(const DoubleReal value) const;

		/// Gets the parameters
		void getParameters(Param& params) const;
		
	protected:
		/// Computes the B-spline fit
		void computeFit_();

		/// Computes the linear extrapolation
    void computeLinear_(const double pos, double& slope, double& offset, 
												double& sd_err);

		/// Vectors for B-spline computation
    gsl_vector *x_, *y_, *w_, *bsplines_, *coeffs_;
		/// Covariance matrix
    gsl_matrix *cov_;
		/// B-spline workspace
    gsl_bspline_workspace *workspace_;
		/// Number of data points and coefficients
    size_t size_, ncoeffs_;
		// First/last breakpoint
    double xmin_, xmax_;
		/// Parameters for linear extrapolation
    double slope_min_, slope_max_, offset_min_, offset_max_;
		/// Fitting errors of linear extrapolation
    double sd_err_left_, sd_err_right_;
	};

} // end of namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_TRANSFORMATIONMODEL_H
