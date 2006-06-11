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
// $Id: SimpleSeeder.h,v 1.17 2006/04/25 12:41:44 ole_st Exp $
// $Author: ole_st $
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

	/** @brief Seeding class that follows the description of Groepl et al (2004)
	  
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
		
	*/ 
  class SimpleSeeder 
    : public BaseSeeder
  {

  public:
		/** @brief Functor that allows to compare the indizes 
				of two peaks by their intensity.

				\internal (cg) Note that std::sort() will copy the comparator object
				over and over again, hence a comparator without data members is
				preferable.
				(ost) But in this context not feasible...
		*/
  	class IntensityLess 
  	{
  	
  		public:
  			/// Construtor
  			IntensityLess(FeaFiTraits* traits)
  				: traits_(traits) 
  			{}
  			
  			/// Overloaded () operator that allows to treat this class as a functor.
  			bool operator() (const Index& x, const Index& y)
				{
    			return traits_->getPeakIntensity(x) < traits_->getPeakIntensity(y);
				}
  			
  		protected:
  			FeaFiTraits* traits_;
  	};
  	
    /// standard constructor
    SimpleSeeder();

    /// destructor 
    virtual ~SimpleSeeder();

    /// return next seed 
    Index nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new SimpleSeeder();
    }

    static const String getName()
    {
      return "SimpleSeeder";
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
  	
  	/// counts the number of seeds that we returned so far
  	unsigned int nr_seeds_;
  	
  	/// the assumed noise threshold as a percentage of the fifth largest peak
  	double noise_threshold_;

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
