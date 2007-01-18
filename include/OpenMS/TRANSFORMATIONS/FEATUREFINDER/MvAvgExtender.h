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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MVAVGEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MVAVGEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/DATASTRUCTURES/RunningAveragePosition.h>

#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/KernelTraits.h>

#include <OpenMS/MATH/MISC/LinearInterpolation.h>

#include <queue>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

namespace OpenMS
{

	/**
	  @brief Implements the extension phase of the FeatureFinder 
	
			@ingroup FeatureFinder
		
	 */
  class MvAvgExtender 
    : public BaseExtender
  {

  public:
  
		/// Intensity of a data point
  	typedef FeaFiTraits::IntensityType IntensityType;
		/// Coordinates of a point (m/z and rt)
  	typedef FeaFiTraits::CoordinateType CoordinateType;
		/// Priority of a point (see below)
  	typedef KernelTraits::ProbabilityType ProbabilityType;
		/// Point type (raw, picked etc)
		typedef FeaFiTraits::PeakType PeakType;
		/// Position of a point
		typedef PeakType::PositionType PositionType;
  	
  	enum DimensionId
			{ 
				RT = DimensionDescription < LCMS_Tag >::RT,
				MZ = DimensionDescription < LCMS_Tag >::MZ
			};
  
  	/// standard constructor
    MvAvgExtender();

    /// destructor
    virtual ~MvAvgExtender();

    /// return next seed
    const IndexSet& extend(const IndexSet& seed_region);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new MvAvgExtender();
    }

		/// returns the name of this module 
    static const String getName()
    {
      return "MvAvgExtender";
    }
    
    /**
     *  @brief A helper structure to sort indizes by their priority.
     * 
     *  This structure is used to keep track of the boundary of a 
     *  feature. After a peak is found during the extension phase,
     *  we compute its priority (which is dependant on its distance from
     *  the point that was the last to be extracted from the boundary
     *  and its intensity). If this priority is large enough, we include
     *  the point into the boundary. The boundary (which is implemented
     *  as mutable priority queue) sorts the peaks by this priority.
     * 
     * */ 	
  	struct IndexWithPriority
  	{
  		IndexWithPriority(UnsignedInt i, double p) : index(i), priority(p) {}
  		
  		UnsignedInt index;
  		ProbabilityType priority;
  		
			/// @brief Compares two indizes by priority.		
  		struct PriorityLess
  		{  			 			  			
  			
  			inline bool operator() (const IndexWithPriority& x, const IndexWithPriority& y) const
				{
    			return x.priority < y.priority;
				}
		
			};
  			
  	};
            
  protected:
  	 	
  	/// Checks if the current peak is too far from the centroid
  	bool isTooFarFromCentroid_(UnsignedInt current_index);
   	
   	/// Extends the seed into positive m/z direction   	
  	void moveMzUp_(UnsignedInt current_peak);
  	
  	/// Extends the seed into negative m/z direction 
  	void moveMzDown_(UnsignedInt current_peak);
  	
  	/// Extension into positive rt dimension 
  	void moveRtUp_(UnsignedInt current_peak);
  	
  	/// Extends the seed into negative retention time direction 
  	void moveRtDown_(UnsignedInt current_peak);
  	
  	/// Computes the priority of a peak as function of intensity and distance from seed. 
  	ProbabilityType computePeakPriority_(PeakType & peak);
  	
  	/// Checks the neighbours of the current for insertion into the boundary.
  	void checkNeighbour_(UnsignedInt current_peak);
  	
  	/// This flag indicates whether the first seed has already been processed. 
  	bool first_seed_seen_;
  	
  	/// Data points with intensity below this threshold are not considered in the extension phase. 
  	IntensityType intensity_threshold_;
  	
		/// Factor for minimum seed intensity 
  	IntensityType intensity_factor_;
  	  	
  	/// keeps an running average of the peak coordinates weighted by the intensities 
  	RunningAveragePosition< PositionType > running_avg_;
  	
  	/// Keeps track of peaks already included in the boundary (value is priority of peak) 
  	std::map<UnsignedInt, ProbabilityType> priorities_; 
  	
  	/// Position of last peak extracted from the boundary (used to compute the priority of neighbouring peaks)
  	PositionType last_pos_extracted_;
		  	
  	/// Represents the boundary of a feature 
  	std::priority_queue< IndexWithPriority, std::vector < IndexWithPriority > , IndexWithPriority::PriorityLess > boundary_;    
  	
		/// Score distribution in retention time     			
  	Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_rt_;

		/// Score distribution in m/z
		Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_mz_;
		
		/// Sum of the intensities collected so far
		IntensityType intensity_sum_;
		
		/// Mininum percentage of the already collected intensity that a new pointed has to contribute
		IntensityType min_intensity_contribution_;
		
		/// Maximum distance to seed in positive m/z
		CoordinateType dist_mz_up_; 
		/// Maximum distance to seed in negative m/z
		CoordinateType dist_mz_down_; 
		/// Maximum distance to seed in positive retention time
		CoordinateType dist_rt_up_; 
		/// Maximum distance to seed in negative retention time
		CoordinateType dist_rt_down_;   			
		
		/// Minium priority for points in the feature region (priority is function of intensity and distance to seed)
		ProbabilityType priority_threshold_;
		
		PositionType seed_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MVAVGEXTENDER_H
