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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FASTEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FASTEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/DATASTRUCTURES/RunningAveragePosition.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ComparatorUtils.h> 
#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

namespace OpenMS
{

	/**
		 @brief Experimental Extender ! Do not use this one !
	  
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
  class FastExtender 
    : public BaseExtender
  {

	 public:
  
  	typedef FeaFiTraits::IntensityType IntensityType;
  	typedef FeaFiTraits::CoordinateType CoordinateType;
  	typedef KernelTraits::ProbabilityType ProbabilityType;
  	
  	enum DimensionId
			{ 
				RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
				MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
			};
  
  	/// standard constructor
    FastExtender();

    /// destructor
    virtual ~FastExtender();

    /// return next seed
    const IndexSet& extend(const UnsignedInt seed);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new FastExtender();
    }

		/// returns the name of this module 
    static const String getName()
    {
      return "FastExtender";
    }
    
    /**
      @brief A helper structure to sort indizes by their priority.
     
      This structure is used to keep track of the boundary of a feature. After
      a peak is found during the extension phase, we compute its priority
      (which is dependant on its distance from the point that was the last to
      be extracted from the boundary and its intensity). If this priority is
      large enough, we include the point into the boundary. The boundary
      (which is implemented as mutable priority queue) sorts the peaks by this
      priority.
			
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
  	 	
  	/// Checks whether the current peak is too far from the centroid
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
  	double computePeakPriority_(UnsignedInt current_peak);
  	
  	/// Checks the neighbours of the current for insertion into the boundary.
  	void checkNeighbour_(UnsignedInt current_peak);
  	
  	/// This flag indicates whether the first seed has already been processed. 
  	bool first_seed_seen_;
    
		/// Tolerance of the sum of intensities in the feature regions	
  	double intensity_factor_;
  	  	
  	/// keeps an running average of the peak coordinates weighted by the intensities 
  	RunningAveragePosition< DPosition<2> > running_avg_;
  	
  	/// Keeps track of peaks already included in the boundary (value is priority of peak) 
  	std::map<UnsignedInt, double> priorities_; 
  	  	
  	/// The last peak that was extracted from the boundary (used to compute the priority of neighbouring peaks)
  	UnsignedInt last_extracted_;
  	
  	/// Represents the boundary of a feature 
  	std::priority_queue<IndexWithPriority, 
  	            				std::vector<IndexWithPriority>, 
  	           					IndexWithPriority::PriorityLess> boundary_;    
  	
		/// Interpolates the priority of a peak in rt	             			
  	Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_rt_;

		/// Interpolates the priority of a peak in mz
		Math::LinearInterpolation < CoordinateType, ProbabilityType > score_distribution_mz_;
		
		/// Counts the number of peak we encountered so far
		UnsignedInt nr_peaks_seen_;
	
		/// The sum of the peaks in the feature region
		IntensityType intensity_sum_;
  	
		/// A moving average of the intensities collected so far
  	std::vector<IntensityType> moving_avg_;
  	
		/// The last moving average we computed 
  	IntensityType last_avg_;
  	
		/// The tolerance of the average intensities
  	IntensityType intensity_avg_tol_;  	
  	
  	float dist_mz_up_; 
		float dist_mz_down_; 
		float dist_rt_up_; 
		float dist_rt_down_;   	
	
		/// Peaks with intensity below this baseline are ignored		
		float extension_baseline_;
	
		/// Initializes the extender class.
		void init_();
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H

