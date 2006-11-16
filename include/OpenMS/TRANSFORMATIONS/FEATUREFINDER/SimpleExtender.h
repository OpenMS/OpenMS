// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/DATASTRUCTURES/RunningAveragePosition.h>

#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ComparatorUtils.h> 
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
					<td>Scale for the interpolation of the rt distribution (default 2.0)</td></tr>
			<tr><td>tolerance_mz</td>
					<td>Scale for the interpolation of the mz distribution (default 0.5)</td></tr>
			<tr><td>dist_mz_up</td>
					<td>Maximum distance in positive mz direction for data points in the feature region (default 6.0)</td></tr>
			<tr><td>dist_mz_down</td>
					<td>Maximum distance in negative mz direction for data points in the feature region (default 2.0)</td></tr>
			<tr><td>dist_rt_up</td>
					<td>Maximum distance in postive rt direction for data points in the feature region (default 5.0)</td></tr>
			<tr><td>dist_rt_down</td>
					<td>Maximum distance in negative rt direction for data points in the feature region (default 5.0)</td></tr>
			<tr><td>priority_thr</td>
					<td>Minimum priority for data points to be included into the boundary of the feature (default 0.01)</td></tr>
			<tr><td>intensity_factor</td>
					<td>Influences for intensity (ion count) threshold in the feature extension. We include only raw data
					points into this region if their intensity is larger than [intensity_factor * (intensity of the 5th largest peak in the region)].
					We use the 5th largest peak for robustness reasons (default 0.03) .</td></tr>
			</table>
			
			@ingroup FeatureFinder
		
	 */
  class SimpleExtender 
    : public BaseExtender
  {

  public:
  
  	typedef FeaFiTraits::IntensityType IntensityType;
  	typedef FeaFiTraits::CoordinateType CoordinateType;
  	typedef KernelTraits::ProbabilityType ProbabilityType;
		typedef FeaFiTraits::PeakType PeakType;
  	
  	enum DimensionId
			{ 
				RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
				MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
			};
  
  	/// standard constructor
    SimpleExtender();

    /// destructor
    virtual ~SimpleExtender();

    /// return next seed
    const IndexSet& extend(const UnsignedInt seed);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new SimpleExtender();
    }

		/// returns the name of this module 
    static const String getName()
    {
      return "SimpleExtender";
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
  		double priority;
  		  		
  		struct PriorityLess
  		{  			 			  			
  			/// Compare two indizes by priority.
  			inline bool operator() (const IndexWithPriority& x, const IndexWithPriority& y) const
				{
    			return x.priority < y.priority;
				}
		
			};
  			
  	};
  	
  	
  	/**
  	 	 @brief This structure simulates a boost::property_map and simply
  	   returns the index for each IndexWithPriority structure.
  	  
  	   Only needed for the priority queue implementing the feature
  	   boundary.
  	 **/
  	struct IndexMap
  	{
  		typedef int value_type;
			typedef IndexWithPriority key_type;
					
  		IndexMap() {}
  		
  		inline int operator[] (const IndexWithPriority& iwp)  const throw()
  		{
  			return iwp.index;
  		}
  		
  		friend int get ( IndexMap const & /*imap*/, IndexWithPriority const & key ) throw()
			{ 
				return key.index; 
			}
  		  		
  	};
          
  protected:
  	 	
  	/// Checks if the current peak is too far from the centroid
  	bool isTooFarFromCentroid_(UnsignedInt& current_index);
   	
   	/// Extends the seed into positive m/z direction   	
  	void moveMzUp_(UnsignedInt current_peak);
  	
  	/// Extends the seed into negative m/z direction 
  	void moveMzDown_(UnsignedInt current_peak);
  	
  	/// Extension into positive rt dimension 
  	void moveRtUp_(UnsignedInt current_peak);
  	
  	/// Extends the seed into negative retention time direction 
  	void moveRtDown_(UnsignedInt current_peak);
  	
  	/// Computes the priority of a peak as function of intensity and distance from seed. 
  	double computePeakPriority_(const PeakType& peak);
  	
  	/// Checks the neighbours of the current for insertion into the boundary.
  	void checkNeighbour_(UnsignedInt& peak_index);
  	
  	/// This flag indicates whether the first seed has already been processed. 
  	bool first_seed_seen_;
  	
  	/// Data points with intensity below this threshold are not considered in the extension phase. 
  	IntensityType intensity_threshold_;
  	
		/// Factor for minimum seed intensity 
  	IntensityType intensity_factor_;
  	  	
  	/// keeps an running average of the peak coordinates weighted by the intensities 
  	RunningAveragePosition< DPosition<2> > running_avg_;
  	
  	/// Keeps track of peaks already included in the boundary (value is priority of peak) 
  	std::map<UnsignedInt, double> priorities_; 
  	
  	/// Position of last peak extracted from the boundary (used to compute the priority of neighbouring peaks)
  	DPosition<2> last_pos_extracted_;
		  	
  	/// Represents the boundary of a feature 
  	std::priority_queue< IndexWithPriority, std::vector < IndexWithPriority > , IndexWithPriority::PriorityLess > boundary_;    
  	           					
  	/*MutablePriorityQueue<IndexWithPriority, 
  	                    std::vector<IndexWithPriority>, 
  	                    IndexWithPriority::PriorityLess,
  	                    IndexMap > boundary_; */ 	            
  	
		/// Score distribution in m/z dimension             			
  	Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_rt_;

		/// Score distribution in time dimension      
		Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_mz_;
		
		/// Number of peaks encountered so far
		UnsignedInt nr_peaks_seen_;
		
		/// Maximum distance to seed in positive m/z
		CoordinateType dist_mz_up_; 
		/// Maximum distance to seed in negative m/z
		CoordinateType dist_mz_down_; 
		/// Maximum distance to seed in positive retention time
		CoordinateType dist_rt_up_; 
		/// Maximum distance to seed in negative retention time
		CoordinateType dist_rt_down_;   	
		
		/// Minium priority for points in the feature region (priority is function of intensity and distance to seed)
		float priority_threshold_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H
