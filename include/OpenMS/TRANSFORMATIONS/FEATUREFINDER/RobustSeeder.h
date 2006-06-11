// -*- C++: make; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: RobustSeeder.h,v 1.5 2006/04/25 12:41:44 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ROBUSTSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ROBUSTSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include<algorithm>
#include<vector>
#include<iostream>
#include<math.h>

namespace OpenMS
{
	/** @brief Seeding module for the feature finding algorithmus.
	 
	 		This seeder simply collects the n largest peaks and computes
	 		the average. A fixed percentage of this average is used as
	 		threshold for the seed intensity. We believe that this
	 		procedure is more robust than simply chosing a percentage
	 		of the fifth largest peak as threshold.
	 		
	 		<table>
			<tr><td>min_intensity</td>
					<td>absolute value for the minimum intensity of a seed
					    if not set, a fixed percentage (see below) of the average intensity of
					    the 300 largest peaks is taken (no default).</td></tr>
			<tr><td>intensity_perc</td>
					<td>percentage (see below) of the average intensity of
					    the 300 largest peaks as minimum seed intensity (default: 0.3 = 30 %).</td></tr>
			</table>
			
			@ingroup FeatureFinder
		
	*/
  class RobustSeeder 
    : public BaseSeeder
  {

  public:
  	
  	typedef FeaFiTraits::IntensityType IntensityType;
  
		/** @brief Functor that allows to compare the indizes 
				of two peaks by their intensity.
				
		*/
  	class IntensityLess 
  	{
  	
  		public:
  			/** @brief Construtor **/
  			IntensityLess(FeaFiTraits* traits) 
  				: traits_(traits)
  			{}
  			
  			/** @brief Overloaded () operator that allows to
  			           treat this class as a functor. */
  			bool operator() (const Index& x, const Index& y)
				{
    			return traits_->getPeakIntensity(x) < traits_->getPeakIntensity(y);
				}
  			
  		protected:
  			FeaFiTraits* traits_;
  	};
  	
    /// standard constructor
    RobustSeeder();

    /// destructor
    virtual ~RobustSeeder();

    /// return next seed 
    Index nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new RobustSeeder();
    }

    static const String getName()
    {
      return "RobustSeeder";
    }

  protected:
    /// sort the indizes according to peak intensity 
  	void sort_();
  	
  	/// contains the indizes 
  	std::vector<UnsignedInt> indizes_;
  	
  	/// Indicates whether the vector of indizes is sorted
  	bool is_initialised_;
  	
  	/// Points to the next peak in the peak vector 
  	std::vector<UnsignedInt>::const_iterator current_peak_;
  	  	
  	/// the assumed noise threshold as a percentage of the fifth largest peak 
   	FeaFiTraits::IntensityType noise_threshold_;
   	
   	/// Number of seeds encountered so far
   	UnsignedInt nr_seeds_;  	
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_ROBUSTSEEDER_H
