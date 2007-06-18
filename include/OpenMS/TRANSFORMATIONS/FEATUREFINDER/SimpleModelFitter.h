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
