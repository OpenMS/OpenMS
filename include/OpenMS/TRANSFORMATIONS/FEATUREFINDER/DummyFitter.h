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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYFITTER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

namespace OpenMS
{
	
	/**
		@brief Dummy fitting. Use this module if you don't need any refinement of your features by model fitting.
		
		This module does not do anything besides constructing the feature from
		the region.
		
		<table>
			<tr>
				<td>min_num_peaks:final</td>
				<td>minimum number of peaks in the feature region.</td>
			</tr>
			<tr>
				<td>min_num_peaks:extended</td>
				<td>minimum number of peaks in the region after extension.</td>
			</tr>
			<tr>
				<td>use_fwhm_intensity</td>
				<td>binary value (0/1). Use fwhm (full witdth at half maximum) for quantification or not..</td>
			</tr>
		</table> 
    
    @ref DummyFitter_Parameters are explained on a separate page.
		 
		@ingroup FeatureFinder
		
  */
  class DummyFitter
    : public BaseModelFitter
  {
  	public:
						
			///
			typedef Feature::CoordinateType FeatureCoordinateType;
			///
			typedef Feature::PositionType PositionType2D;
			///
			typedef FeaFiTraits::CoordinateType CoordinateType;
			///
			typedef FeaFiTraits::IntensityType IntensityType;
			
	    /// Default constructor
	    DummyFitter();
	
	    /// destructor
	    virtual ~DummyFitter();

	    /// Copy constructor
	    DummyFitter(const DummyFitter& rhs);
	    
	    /// Assignment operator
	    DummyFitter& operator= (const DummyFitter& rhs);
	
	    /// return next seed
	    Feature fit(const ChargedIndexSet& range) throw (UnableToFit);
	
	    static BaseModelFitter* create()
	    {
	      return new DummyFitter();
	    }
	
	    static const String getProductName()
	    {
	      return "DummyFitter";
	    }
			
			protected:
				UInt counter_;
	
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYFITTER_H
