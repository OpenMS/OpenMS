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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEMODELFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEMODELFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/MATH/STATISTICS/AsymmetricStatistics.h>

namespace OpenMS
{
	class BaseQuality;

	/**
		@brief Matches a peptide template (based on averagines) to signals extracted from a LC-MS map.
		
		The template is two-dimensional (consists of a part for m/z and rt). The m/z domain is
		based on averagines and also contains a parameter for the mass (in-)accuracy of the
		mass analyser ("isotope_model:stdev"). To achieve a higher performance, you should set the
		range of tested value to a small value.
			 
		Parameters:
		<table>
		<tr><td></td><td></td><td>tolerance_stdev_bounding_box</td>
		<td>bounding box has range [minimim of data, maximum of data] enlarged
		by tolerance_stdev_bounding_box times the standard deviation of the data</td></tr>
		<tr><td></td><td></td><td>feature_intensity_max</td>
		<td>If this parameter is set to one (default) the peptide abundance is estimated as the
		sum of all peak intensities within the feature region. Otherwise this abundance is
		estimated as the maximum intensity.</td></tr>
		<tr><td></td><td></td><td>intensity_cutoff_factor</td>
		<td>cutoff peaks with a predicted intensity below intensity_cutoff_factor times the maximal intensity of the model</td></tr>
		<tr><td></td><td></td><td>feature_intensity_sum</td>
		<td>estimate abundance of the compound as the sum of all contained peaks</td></tr>
		<tr><td colspan="2">rt</td><td>interpolation_step</td>
		<td>step size in seconds used to interpolate model for rt</td></tr>
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
		 
		@ref SimpleModelFitter_Parameters are explained on a separate page.
		
		@ingroup FeatureFinder
  */
  class SimpleModelFitter
    : public BaseModelFitter
  {
	 public:
	 
	 	enum 
			{
				RT = RawDataPoint2D::RT,
				MZ = RawDataPoint2D::MZ
			};
	  
			/// Ion count
			typedef FeaFiTraits::IntensityType IntensityType;
			/// Single coordinate
			typedef FeaFiTraits::CoordinateType Coordinate;
			///	Single coordinate
			typedef Feature::CoordinateType CoordinateType;
			/// Position in 2D
			typedef Feature::PositionType PositionType2D;
			/// Quality of a feature
			typedef Feature::QualityType QualityType;
			
			enum RtFitting { RTGAUSS=0, BIGAUSS=1};
			
			enum MzFitting { MZGAUSS=0, CHARGE1=1, CHARGE2=2, CHARGE3=3, CHARGE4=4	, CHARGE5=5, CHARGE6=6 };
	
	    /// Default constructor
	    SimpleModelFitter();
	
	    /// destructor
	    virtual ~SimpleModelFitter();

	    /// Copy constructor
	    SimpleModelFitter(const SimpleModelFitter& rhs);
	    
	    /// Assignment operator
	    SimpleModelFitter& operator= (const SimpleModelFitter& rhs);
	
	    /// return next seed
	    Feature fit(const ChargedIndexSet& range) throw (UnableToFit);
	
	    static BaseModelFitter* create()
	    {
	      return new SimpleModelFitter();
	    }
	
	    static const String getProductName()
	    {
	      return "SimpleModelFitter";
	    }
	
		protected:
			virtual void updateMembers_();
			
			/// fit offset by maximizing of quality
			double fitOffset_(	InterpolationModel* model, const IndexSet& set, double stdev1, double stdev2, Coordinate offset_step);
	
			double fit_(	const IndexSet& set, MzFitting mz_fit, RtFitting rt_fit, Coordinate isotope_stdev=0.1);
	
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
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEMODELFITTER_H
