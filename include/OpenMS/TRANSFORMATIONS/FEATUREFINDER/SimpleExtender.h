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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/DATASTRUCTURES/RunningAveragePosition.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>

#include <queue>

namespace OpenMS
{

	/**
	  @brief Implements the extension phase of the FeatureFinder as described by Groepl et al. (2005)
	  
		We want to determine a region around a seed that is
		provided by the seeder. Initially, this region is
		empty. The boundary of this region is implemented
		using a MutablePriorityQueue which contains only
		the seed at the beginning.

		At each step, we choose a data point from the boundary,
		move it into the region and explore the neigbourhood of
		this point in a cross-wise manner (m/z up, m/z down, rt up
    and rt down). During this exploration we compute the priority
		of all encountered points as a function of the distance from
		the extracted point. If this priority exceeds a threshold,
		we insert the corresponding point into the boundary and proceed.

		We stop the extension phase if all peaks contained in the
    boundary have an intensity lower than a threshold or are too
    distant from the centroid of the feature.

		Parameters:
			<table>
			<tr><td>tolerance_rt</td>
					<td>Scale for the interpolation of the rt priority distribution (default 2.0)</td></tr>
			<tr><td>tolerance_mz</td>
					<td>Scale for the interpolation of the mz priority distribution (default 0.5)</td></tr>
			<tr><td>dist_mz_up</td>
					<td>Maximum distance in positive mz direction for data points in the feature region (default 6.0)</td></tr>
			<tr><td>dist_mz_down</td>
					<td>Maximum distance in negative mz direction for data points in the feature region (default 2.0)</td></tr>
			<tr><td>dist_rt_up</td>
					<td>Maximum distance in postive rt direction for data points in the feature region (default 5.0)</td></tr>
			<tr><td>dist_rt_down</td>
					<td>Maximum distance in negative rt direction for data points in the feature region (default 5.0)</td></tr>
			<tr><td>priority_thr</td>
					<td>Minimum priority for data points to be included into the boundary of the feature (default 0.0)
								The priority of a data point is a function of its intensity and its distance to the last point
								included into the feature region. Setting this threshold to zero or a very small value is
								usually a good idea. </td></tr>
			<tr><td>intensity_factor</td>
					<td>Influences for intensity (ion count) threshold in the feature extension. We include only raw data
					points into this region if their intensity is larger than [intensity_factor * (intensity of the seed)].
					(default value is 0.03) .</td></tr>
			</table>
		
		@todo Test on different data types (peak, raw, low/high intensity) (Ole)
		@todo Use priorities or remove them (Ole)
		@todo Try to divide intensity_sum by region_.size() (Ole)
		
		@ingroup FeatureFinder
	*/
  class SimpleExtender 
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
    SimpleExtender();

    /// destructor
    virtual ~SimpleExtender();

    /// Copy constructor
    SimpleExtender(const SimpleExtender& rhs);
    
    /// Assignment operator
    SimpleExtender& operator= (const SimpleExtender& rhs);

    /// return next seed
    const IndexSet& extend(const IndexSet& seed_region);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new SimpleExtender();
    }

		/// returns the name of this module 
    static const String getProductName()
    {
      return "SimpleExtender";
    }
    
    /**
     @brief A helper structure to sort indizes by their priority.
     
     This structure is used to keep track of the boundary of a 
     feature. After a peak is found during the extension phase,
     we compute its priority (which is dependant on its distance from
     the point that was the last to be extracted from the boundary
     and its intensity). If this priority is large enough, we include
     the point into the boundary. The boundary (which is implemented
     as mutable priority queue) sorts the peaks by this priority.
     
    */ 	
  	struct IndexWithPriority
  	{
  		IndexWithPriority(const IDX& i, double p) : index(i), priority(p)
  		{
  		}
  		
  		IDX index;
  		ProbabilityType priority;
  		
			///Compares two indizes by priority.		
  		struct PriorityLess
  		{  			 			  			
  			
  			inline bool operator() (const IndexWithPriority& x, const IndexWithPriority& y) const
				{
    			return x.priority < y.priority;
				}
		
			};
  			
  	};
            
  protected:
  	virtual void updateMembers_();
  	
  	/// Checks if the current peak is too far from the centroid
  	bool isTooFarFromCentroid_(const IDX& current_index);
   	
   	/// Extends the seed into positive m/z direction   	
  	void moveMzUp_(const IDX& current_peak);
  	
  	/// Extends the seed into negative m/z direction 
  	void moveMzDown_(const IDX& current_peak);
  	
  	/// Extension into positive rt dimension 
  	void moveRtUp_(const IDX& current_peak);
  	
  	/// Extends the seed into negative retention time direction 
  	void moveRtDown_(const IDX& current_peak);
  	
  	/// Computes the priority of a peak as function of intensity and distance from seed. 
  	ProbabilityType computePeakPriority_(const IDX& index);
  	
  	/// Checks the neighbours of the current for insertion into the boundary.
  	void checkNeighbour_(const IDX& index);
  	  	
  	/// keeps an running average of the peak coordinates weighted by the intensities 
  	RunningAveragePosition< PositionType2D > running_avg_;
  	
  	/// Keeps track of peaks already included in the boundary (value is priority of peak) 
  	std::map<IDX, ProbabilityType> priorities_; 
  	
  	/// Position of last peak extracted from the boundary (used to compute the priority of neighbouring peaks)
  	PositionType2D last_pos_extracted_;
		  	
  	/// Represents the boundary of a feature 
  	std::priority_queue< IndexWithPriority, std::vector < IndexWithPriority > , IndexWithPriority::PriorityLess > boundary_;    
  	
		/// Score distribution in retention time     			
  	Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_rt_;

		/// Score distribution in m/z
		Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_mz_;
		
		/// Mininum intensity of a boundary point. Calculated from 'intensity_factor' and the seed intensity
		IntensityType intensity_threshold_;
		
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
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H
