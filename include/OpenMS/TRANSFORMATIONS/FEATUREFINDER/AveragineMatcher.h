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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_AVERAGINEMATCHER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_AVERAGINEMATCHER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>

#include <OpenMS/MATH/STATISTICS/AsymmetricStatistics.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

namespace OpenMS
{

	class BaseQuality;

	/**
		 @brief Extended model fitter using gaussian or isotope model in mz and bigauss, lmagauss (bigauss with parameters estimated by levenberg-marquardt) 
		 or emg (exponent. modified gaussian) for the retention time domain.

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
		 
		 @ingroup FeatureFinder
  */


  class AveragineMatcher
    : public BaseModelFitter
  {

	 public:
		///
		typedef IndexSet::const_iterator IndexSetIter;
		///
		typedef FeaFiTraits::CoordinateType Coordinate;
		///
		typedef Feature::CoordinateType CoordinateType;
		/// 
		typedef Feature::QualityType QualityType;
		///	
		typedef Feature::IntensityType IntensityType;
		///	
		typedef Feature::PositionType PositionType2D;
		///
		typedef Feature::ChargeType ChargeType;
		
				
		enum RtFitting{ RTGAUSS=0, LMAGAUSS=1, EMGAUSS=2, BIGAUSS=3, LOGNORMAL=4 };
		enum MzFitting{ MZGAUSS=0, CHARGE1=1, CHARGE2=2, CHARGE3=3, CHARGE4=4	};
		
		

		enum 
			{
				RT = RawDataPoint2D::RT,
				MZ = RawDataPoint2D::MZ
			};

    /// Default constructor
    AveragineMatcher();

    /// Destructor
    virtual ~AveragineMatcher();

   	/// Copy constructor
    AveragineMatcher(const AveragineMatcher& rhs);
    
    /// Assignment operator
    AveragineMatcher& operator= (const AveragineMatcher& rhs);

    /// Return next feature
    Feature fit(const ChargedIndexSet& range) throw (UnableToFit);

    static BaseModelFitter* create()
    {
      return new AveragineMatcher();
    }

    static const String getProductName()
    {
      return "AveragineMatcher";
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
	 
	 /// debug info
	 	void dump_all_(ChargedIndexSet set, UInt sampling_size);
	 
		/// loop
	 QualityType fit_loop_(const ChargedIndexSet& set, Int& first_mz, Int& last_mz, CoordinateType& sampling_size_mz , ProductModel<2>*& final); 

		virtual void updateMembers_();

		QualityType fitOffset_(InterpolationModel* model, const IndexSet& set, double stdev1, double stdev2, Coordinate offset_step);
		
		QualityType fit_mz_(ChargedIndexSet set, UInt samplingsize, MzFitting charge,Coordinate isotope_stdev);

		QualityType fit_(const ChargedIndexSet& set, MzFitting mz_fit, RtFitting rt_fit, Coordinate isotope_stdev, UInt samplingsize);
		
		QualityType compute_mz_corr_(IntensityType& mz_data_sum, IsotopeModel& iso_model, CoordinateType& mz_data_avg);

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
		UInt max_iteration_;

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
		
		/// projection of points onto mz
		Math::LinearInterpolation<CoordinateType,CoordinateType> mz_lin_int_;
		
		/// Averagine template for m/z
		IsotopeModel mz_model_;
		
		/// Exponentially modified Gaussian for retention time
		EmgModel rt_model_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_AVERAGINEMATCHER_H
