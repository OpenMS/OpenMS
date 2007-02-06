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
 		
 		<table>
			<tr>
				<td>min_intensity</td>
				<td>Absolute value for the minimum intensity of a seed. If set to 0, a fixed 
					  percentage of the intensity of the largest peak is taken (see below).</td>
			</tr>
			<tr>
				<td>intensity_perc</td>
				<td>Minimum percentage of the intensity of the largest peak that a seed has to have
				    (used only if min_nitensity is set to 0.</td>
			</tr>
		</table>
		
		@ingroup FeatureFinder
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
					: traits_(traits)
  			{
				}
  			
  			inline bool operator() (const IDX& x, const IDX& y)
				{
    			return traits_->getPeakIntensity(x) < traits_->getPeakIntensity(y);
				}
  			
  		protected:
				FeaFiTraits* traits_;
  	};
  	
    /// Default constructor
    SimpleSeeder();

    /// destructor 
    virtual ~SimpleSeeder();

	  /// Copy constructor
	  SimpleSeeder(const SimpleSeeder& rhs);
	  
	  /// Assignment operator
	  SimpleSeeder& operator= (const SimpleSeeder& rhs);
  
    /// return next seed 
    IndexSet nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new SimpleSeeder();
    }

    static const String getProductName()
    {
      return "SimpleSeeder";
    }

  protected:
  	/// contains the indizes 
  	std::vector<IDX> indizes_;

  	/// Indicates whether the vector of indizes is sorted 
  	bool is_initialized_;
  	
  	/// Points to the next peak in the peak vector 
  	std::vector<IDX>::const_iterator current_peak_;
  	
  	/// counts the number of seeds that we returned so far
  	UnsignedInt nr_seeds_;
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
