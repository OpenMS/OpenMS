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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <queue>
#include <algorithm>
#include <cmath>

namespace OpenMS
{

	/**
	  @brief A straightforward implementation of the extension phase of the feature detection / quantification algorithm in OpenMS.
		
		As one can imagine, this module implements a rather trivial extension.
		
		This module takes a single data point (seed) or a group of points (seeding region) 
		and adds further points until :
		
		 (a) until a given maximum distance from the highest data point in the seeding region is reached or
		 (b) no more points are found in the neighbourhood of the seed that contribute a significant amount 
		 of intensity to the feature.
		 		 
		 <table>
			<tr>
				<td>dist_mz_down / dist_mz_up</td>
				<td>Lower and upper bound for the distance in m/z from the highest point.</td>
			</tr>
			<tr>
				<td>dist_rt_down / dist_rt_up</td>
				<td>Lower and upper bound for the distance in rt from the highest point.</td>
			</tr>
			<tr>
				<td>intensity_perc</td>
				<td>Minimum percentage of the intensity of the largest peak that a seed has to have
				    (used only if min_nitensity is set to 0.</td>
			</tr>
				<td>min_intensity_contribution_</td>
				<td>Minimum percentage of the feature intensity sum collected so far 
				      that a point has to contribute.</td>
			</tr>
		</table>
		
		
		@ingroup FeatureFinder
	*/
  class DummyExtender 
    : public BaseExtender
  {

  public:
  
		/// Intensity of a data point
  	typedef DoubleReal IntensityType;
		/// Coordinates of a point (m/z and rt)
  	typedef DoubleReal CoordinateType;
		/// Priority of a point (see below)
  	typedef DoubleReal ProbabilityType;
		/// Position of a point
		typedef  FeaFiTraits::PositionType2D PositionType2D;
  	/// Default constructor
    DummyExtender();

    /// destructor
    virtual ~DummyExtender();

    /// Copy constructor
    DummyExtender(const DummyExtender& rhs);
    
    /// Assignment operator
    DummyExtender& operator= (const DummyExtender& rhs);

    /// return next seed
    const IndexSet& extend(const IndexSet& seed_region);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new DummyExtender();
    }

		/// returns the name of this module 
    static const String getProductName()
    {
      return "DummyExtender";
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
		/// The name speaks for itself
		IntensityType min_intensity_contribution_;
				
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYEXTENDER_H
