// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Marcel Grunert $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDMODELFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDMODELFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/AsymmetricStatistics.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

namespace OpenMS
{

	class BaseQuality;

	/**
		@brief Extended model fitter using gaussian or isotope model in mz and bigauss, lmagauss (bigauss with Levenberg-Marquardt aproximized parameters) or emg (exponent. modified Gaussian with lma aproximized parameters) in rt.
		
		For the isotope model different charges and deviations are tested.<br>
    
    @ref ExtendedModelFitter_Parameters are explained on a separate page.
		
		@todo Merge ExtendedModelFitter and SimpleModelFitter (Clemens)
		@todo Check use of Enums for RT and m/z fit. They destroy the factory concept! (Clemens)
		
		@ingroup FeatureFinder
  */
  class ExtendedModelFitter
    : public BaseModelFitter
  {

	 public:

		typedef IndexSet::const_iterator IndexSetIter;
		typedef FeaFiTraits::CoordinateType Coordinate;

		typedef Feature::CoordinateType CoordinateType;
		typedef Feature::PositionType PositionType2D;

		enum RtFitting{ RTGAUSS=0, LMAGAUSS=1, EMGAUSS=2, BIGAUSS=3, LOGNORMAL=4 };
		enum MzFitting{ MZGAUSS=0, CHARGE1=1, CHARGE2=2, CHARGE3=3, CHARGE4=4	};

		enum 
			{
				RT = RawDataPoint2D::RT,
				MZ = RawDataPoint2D::MZ
			};

    /// Default constructor
    ExtendedModelFitter();

    /// Destructor
    virtual ~ExtendedModelFitter();

   	/// Copy constructor
    ExtendedModelFitter(const ExtendedModelFitter& rhs);
    
    /// Assignment operator
    ExtendedModelFitter& operator= (const ExtendedModelFitter& rhs);

    /// Return next feature
    Feature fit(const ChargedIndexSet& range) throw (UnableToFit);

    static BaseModelFitter* create()
    {
      return new ExtendedModelFitter();
    }

    static const String getProductName()
    {
      return "ExtendedModelFitter";
    }

	/// create a vector with RT-values & Intensities and compute the parameters (intial values) for the EMG, Gauss and logNormal function
	void setData (const IndexSet& set);

	/// Evaluation of the target function for nonlinear optimization.
	int residual(const gsl_vector* x, void* /* params */, gsl_vector* f);

	/// Compute the Jacobian of the residual, where each row of the matrix corresponds to a point in the data.
	int jacobian(const gsl_vector* x, void* /* params */, gsl_matrix* J);

	/// Driver function for the evaluation of function and jacobian.
	int evaluate(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J);

	/// perform a nonlinear optimization
	void optimize();

	/// get height for the EMG and logNormal model
	CoordinateType getHeight() const;

	/// get width for the EMG and logNormal model
	CoordinateType getWidth() const;

	/// get symmetry for the EMG and logNormal model
	CoordinateType getSymmetry() const;

	/// get retention time for the EMG and logNormal model
	CoordinateType getRT() const;

	/// get standard deviation for the Gauss
	CoordinateType getStandardDeviation() const;

	/// get expected value for the Gauss
	CoordinateType getExpectedValue() const;

	/// get scale factor for the Gauss
	CoordinateType getScaleFactor() const;

	/// get GSL status
	std::string getGSLStatus() const;

	 protected:

		virtual void updateMembers_();

		/// fit offset by maximizing of quality
		double fitOffset_(InterpolationModel* model, const IndexSet& set, double stdev1, double stdev2, Coordinate offset_step);

		double fit_(const IndexSet& set, MzFitting mz_fit, RtFitting rt_fit, Coordinate isotope_stdev=0.1);

		BaseQuality* quality_;
		ProductModel<2> model2D_;
		Math::BasicStatistics<> mz_stat_;
		Math::AsymmetricStatistics<> rt_stat_;
		double stdev_mz_;
		double stdev_rt1_;
		double stdev_rt2_;
		PositionType2D min_;
		PositionType2D max_;

		/// counts features (used for debug output only)
		UInt counter_;
			
		/// interpolation step size (in m/z)
		Coordinate interpolation_step_mz_;
		/// interpolation step size (in retention time)
		Coordinate interpolation_step_rt_;
			
		/// first stdev
		float iso_stdev_first_;
		/// last stdev
		float iso_stdev_last_;
		/// step size
		float iso_stdev_stepsize_;
			
		/// first mz model (0=Gaussian, 1....n = charge )
		Int first_mz_model_;			
		/// last mz model
		Int last_mz_model_;

		/// Maximum number of iterations
		unsigned int max_iteration_;

		/// parameter of log normal function:
		/// r is the ratio between h and the height at which w and s are computed
		double r_;

		/// parameter of emg and log normal function:height
		double height_;
		/// parameter of emg and log normal function: width
		double width_;
		/// parameter of emg and log normal function: symmetry
		double symmetry_;
		/// parameter of emg and log normal function: retention time
		double retention_;
		/// parameter indicates symmetric peaks
		bool symmetric_;
		/// gsl status
		std::string gsl_status_;
		/// function for fitting
		std::string profile_;

		/** Test for the convergence of the sequence by comparing the last iteration step dx with the absolute error epsabs and relative error epsrel to the current position x */
		/// absolute error
		double eps_abs_;
		/// relative error
		double eps_rel_;

		/// parameter of gauss function: standard deviation
		double standard_deviation_;
		/// parameter of gauss function: scale factor
		double scale_factor_;
		/// parameter of gauss function: expected value
		double expected_value_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDMODELFITTER_H
