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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

namespace OpenMS
{
	
	/**
		 @brief Dummy fitting.
		 
		 This module does not do anything (besides constructing the feature).
		 
		 @ingroup FeatureFinder
		
  */
  class PeakFitter
    : public BaseModelFitter
  {
  	public:
			
			enum DimensionId
				{
					RT = DimensionDescription < LCMS_Tag >::RT,
					MZ = DimensionDescription < LCMS_Tag >::MZ
				};
			
			///
			typedef DFeature<2>::CoordinateType FeatureCoordinateType;
			///
			typedef DFeature<2>::PositionType PositionType2D;
			///
			typedef FeaFiTraits::CoordinateType CoordinateType;
			///
			typedef FeaFiTraits::IntensityType IntensityType;
			
	    /// Default constructor
	    PeakFitter();
	
	    /// destructor
	    virtual ~PeakFitter();

	    /// Copy constructor
	    PeakFitter(const PeakFitter& rhs);
	    
	    /// Assignment operator
	    PeakFitter& operator= (const PeakFitter& rhs);
	
	    /// return next seed
	    DFeature<2> fit(const IndexSet& range) throw (UnableToFit);
	
	    static BaseModelFitter* create()
	    {
	      return new PeakFitter();
	    }
	
	    static const String getProductName()
	    {
	      return "PeakFitter";
	    }
			
			protected:
				UnsignedInt counter_;
	
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKFITTER_H
