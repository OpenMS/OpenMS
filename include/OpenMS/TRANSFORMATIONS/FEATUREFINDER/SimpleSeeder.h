// -*- C++: make; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <vector>
#include <iostream>

namespace OpenMS
{

	/** 
		@brief Seeding class that follows the description of Groepl et al (2004)
	  
 		The intensity threshold for seeding is taken to be a fixed percentage
 		of the 5th largest peak.	  
 		
 		<table>
		<tr><td>min_intensity</td>
				<td>absolute value for the minimum intensity of a seed
				    if not set, a fixed percentage (see below) of the intensity of
				    the fifth largest peak is taken (no default).</td></tr>
		<tr><td>intensity_perc</td>
				<td>percentage (see below) of the intensity of
				    the fifth largest peak is taken (default: 0.3 = 30%).</td></tr>
		</table>
		
		@ingroup FeatureFinder
		
		@todo Write test with more than one scan and test not only the intensity (Ole)
	*/ 
  class SimpleSeeder 
    : public BaseSeeder
  {
		typedef FeaFiTraits::IntensityType IntensityType;
		typedef FeaFiTraits::MapType MapType;
				
		enum DimensionId
    {
        RT = DimensionDescription < LCMS_Tag >::RT,
        MZ = DimensionDescription < LCMS_Tag >::MZ
    };

  public:
		///Functor that allows to compare the indizes of two peaks by their intensity.
  	class IntensityLess 
  	{			
  		public:
  			IntensityLess(FeaFiTraits* traits)
					: map_(&(traits->getData()))
  			{
				}
  			
  			inline bool operator() (const IDX& x, const IDX& y)
				{
    			return map_->operator[](x.first)[x.second].getIntensity() < map_->operator[](y.first)[y.second].getIntensity();
				}
  			
  		protected:
				FeaFiTraits::MapType* map_;
  	};
  	
    /// Default constructor
    SimpleSeeder();

    /// destructor 
    virtual ~SimpleSeeder();

    /// return next seed 
    IndexSet nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new SimpleSeeder();
    }

    static const String getName()
    {
      return "SimpleSeeder";
    }

  protected:
  	/// contains the indizes 
  	std::vector<IDX> indizes_;

  	/// Indicates whether the vector of indizes is sorted 
  	bool is_initialised_;
  	
  	/// Points to the next peak in the peak vector 
  	std::vector<IDX>::const_iterator current_peak_;
  	
  	/// counts the number of seeds that we returned so far
  	UnsignedInt nr_seeds_;
  	
  	/// the assumed noise threshold as a percentage of the fifth largest peak
  	IntensityType noise_threshold_;

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
