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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <queue>
#include <algorithm>
#include <cmath>

namespace OpenMS
{

	/**
	  @brief Implements the extension phase.
				
		Seeds are extended 
		 (a) until a given maximum distance from seed is reached,
		 (b) as long as they contribute significantly to the feature intensity.
		
		This class we developed for a third-party project and we should discuss whether we want
		it to become part of OpenMS.		
		
		@ingroup FeatureFinder
	*/
  class PeakExtender 
    : public BaseExtender
  {

  public:
  
		/// Intensity of a data point
  	typedef FeaFiTraits::IntensityType IntensityType;
		/// Coordinates of a point (m/z and rt)
  	typedef FeaFiTraits::CoordinateType CoordinateType;
		/// Priority of a point (see below)
  	typedef KernelTraits::ProbabilityType ProbabilityType;
		/// Position of a point
		typedef  FeaFiTraits::PositionType2D PositionType2D;
  	
  	enum DimensionId
			{ 
				RT = DimensionDescription < LCMS_Tag >::RT,
				MZ = DimensionDescription < LCMS_Tag >::MZ
			};
  
  	/// Default constructor
    PeakExtender();

    /// destructor
    virtual ~PeakExtender();

    /// Copy constructor
    PeakExtender(const PeakExtender& rhs);
    
    /// Assignment operator
    PeakExtender& operator= (const PeakExtender& rhs);

    /// return next seed
    const IndexSet& extend(const IndexSet& seed_region);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new PeakExtender();
    }

		/// returns the name of this module 
    static const String getProductName()
    {
      return "PeakExtender";
    }
                  
  protected:
		///
  	virtual void updateMembers_();
  	
  	/// Checks if the current peak is too far from the centroid
  	bool isTooFarFromSeed_(const IDX& current_index);
   	
   	/// Extends the seed into positive m/z direction   	
  	void moveMzUp_(const IDX& current_peak);
  	
  	/// Extends the seed into negative m/z direction 
  	void moveMzDown_(const IDX& current_peak);
  	
  	/// Extension into positive rt dimension 
  	void moveRtUp_(const IDX& current_peak);
  	
  	/// Extends the seed into negative retention time direction 
  	void moveRtDown_(const IDX& current_peak);
     		
		/// Represents the boundary of a feature 
  	std::priority_queue< IDX > boundary_;    
  	
  	/// Seed position (seed is in this case the point with the highest intensity in the seeding region)
  	PositionType2D seed_pos_;
		  	
		/// Maximum distance to seed in positive m/z
		CoordinateType dist_mz_up_; 
		/// Maximum distance to seed in negative m/z
		CoordinateType dist_mz_down_; 
		/// Maximum distance to seed in positive retention time
		CoordinateType dist_rt_up_; 
		/// Maximum distance to seed in negative retention time
		CoordinateType dist_rt_down_;   			
		
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKEXTENDER_H
