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
		 @brief Extended model fitter using gaussian or isotope model in mz and bigauss, lmagauss (bigauss with levenberg-marquardt aproximized parameters) or emg (exponent. modified gaussian with lma aproximized parameters) in rt.

		 For the isotope model different charges and deviations are tested.<br>
		 Parameters:
		 <table>
		 <tr><td></td><td></td><td>tolerance_stdev_bounding_box</td>
		 <td>bounding box has range [minimim of data, maximum of data] enlarged
		 by tolerance_stdev_bounding_box times the standard deviation of the data</td></tr>
		 <tr><td></td><td></td><td>feature_intensity_sum</td>
		 <td>estimate abundance of the compound as the sum of all contained peaks</td></tr>
		 <tr><td></td><td></td><td>intensity_cutoff_factor</td>
		 <td>cutoff peaks with a predicted intensity below intensity_cutoff_factor times the maximal intensity of the model</td></tr>
		 <tr><td rowspan="5" colspan="2">rt</td><td>interpolation_step</td>
		 <td>step size in seconds used to interpolate model for rt</td></tr>
		 <tr><td>max_iteration</td>
		 <td>maximum number of iteration used by the Levenberg-Marquadt algorithms</td></tr>
		 <tr><td>deltaAbsError</td>
		 <td>absolute error used by the Levenberg-Marquardt algorithms</td></tr>
		 <tr><td>deltaRelError</td>
		 <td>relative error used by the Levenberg-Marquardt algorithms</td></tr>
		 <tr><td>profile</td>
		 <td>type of model to try out in rt, possible are GaussModel,							LmaGaussModel, EmgModel or LogNormalModel</td></tr>
		 <tr><td rowspan="2">mz</td><td></td><td>interpolation_step</td>
		 <td>step size in Thomson used to interpolate model for mz</td></tr>
		 <tr><td>model_type</td><td>first, last</td>
		 <td>first (last) type of model to try out in mz,
		 0 = GaussModel,<br>
		 1 = IsotopeModel with charge +1, ..., <br>
		 n = IsotopeModel with charge +n</td></tr>
		 <tr><td rowspan="2" colspan="2">quality</td><td>type</td>
		 <td>name of class derived from BaseModel, measurement for quality of fit</td></tr>
		 <tr><td>minimum</td>
		 <td>minimum quality of feature, if smaller feature will be discarded</td></tr>
		 <tr><td rowspan="2" colspan="2">min_num_peaks</td><td>extended</td>
		 <td>minimum number of peaks gathered by the BaseExtender.
		 If smaller, feature will be discarded </td></tr>
		 <tr><td>final</td>
		 <td>minimum number of peaks left after cutoff.
		 If smaller, feature will be discarded.</td></tr>
		 <tr><td rowspan="5">isotope_model</td><td>stdev</td><td>first, last, step</td>
		 <td>testing isotope standard deviations in range [stdev_first_mz,stdev_last_mz]
		 in steps of size stdev_step_mz.
		 Used to account for different data resolutions</td></tr>
		 <tr><td>avergines</td><td>C, H, N, O, S</td>
		 <td>averagines are used to approximate the number of atoms of a given element
		 (C,H,N,O,S) given a mass</td></tr>
		 <tr><td rowspan="3">isotope</td><td>trim_right_cutoff</td>
		 <td>use only isotopes with abundancies above this cutoff</td></tr>
		 <tr><td>maximum</td>
		 <td>maximum number of isotopes being used for the IsotopeModel</td></tr>
		 <tr><td>distance</td>
		 <td>distance between two isotopes of charge +1</td></tr>
		 </table>
		 
		 @todo use DefaultParamHandler (Clemens, Marcel)
  */


  class ExtendedModelFitter
    : public BaseModelFitter
  {

	 public:

		typedef IndexSet::const_iterator IndexSetIter;
		typedef FeaFiTraits::CoordinateType Coordinate;

		enum DimensionId
			{
				RT = DimensionDescription < LCMS_Tag >::RT,
				MZ = DimensionDescription < LCMS_Tag >::MZ
			};

		typedef DFeature<2>::CoordinateType CoordinateType;
		typedef DFeature<2>::PositionType PositionType2D;

		enum RtFitting{ RTGAUSS=0, LMAGAUSS=1, EMGAUSS=2, BIGAUSS=3, LOGNORMAL=4 };
		enum MzFitting{ MZGAUSS=0, CHARGE1=1, CHARGE2=2, CHARGE3=3, CHARGE4=4	};


    /// Default constructor
    ExtendedModelFitter();

    /// destructor
    virtual ~ExtendedModelFitter();

   	/// Copy constructor
    ExtendedModelFitter(const ExtendedModelFitter& rhs);
    
    /// Assignment operator
    ExtendedModelFitter& operator= (const ExtendedModelFitter& rhs);

    /// return next seed
    DFeature<2> fit(const IndexSet& range) throw (UnableToFit);

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

	 protected:

		virtual void updateMembers_();

		/// fit offset by maximizing of quality
		double fitOffset_(InterpolationModel<>* model, const IndexSet& set, const double stdev1, const double stdev2, const Coordinate offset_step);

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

		unsigned int counter_;

		/// Maximum number of iterations
		unsigned int max_iteration_;

		/// parameter of log normal function
		/// r is the ratio between h and the height at which w and s are computed
		double r_;

		/// parameter of emg and log normal function
		double height_;
		double width_;
		double symmetry_;
		double retention_;
		double deviation_;

		/// the name of the function
		std::string profile_;

		/** 
			Test for the convergence of the sequence by comparing the last iteration step dx
    	with the absolute error epsabs and relative error epsrel to the current position x 
    */
		double eps_abs_;
		double eps_rel_;

		/// parameter of gauss function
		double standard_deviation_;
		double scale_factor_;
		double expected_value_;

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDMODELFITTER_H
