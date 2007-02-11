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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PICKEDPEAKSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PICKEDPEAKSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>

namespace OpenMS
{
	/** 
		@brief Seeding module which tries to find seeds by looking at 

		This extender module sweeps through the scans and classifies cluster 
 		of peaks as candidate peptides if the distance between successive peaks 
 		is 1 Da (charge 1) , 0.5 Da (charge 2) or 0.3 Da (charge 3) etc.
 
 		@note Experiments have shown that this extender produces a lot of false positive hits. It would be
 		better to check if the relaitive intensities between the peaks is similar to an isotopic pattern.
  	
  	@todo Write test with more than one scan and test not only the intensity (Ole)
		@todo Document parameters (Ole)
		  	
		@ingroup FeatureFinder
	*/
	class PickedPeakSeeder
	  : public BaseSeeder
	{
		public:
		  typedef FeaFiTraits::IntensityType IntensityType;
		  typedef FeaFiTraits::CoordinateType CoordinateType;
		  typedef KernelTraits::ProbabilityType ProbabilityType;
			typedef std::multimap<CoordinateType,IsotopeCluster> TableType;
		
		  enum DimensionId
		  {
	      RT = DimensionDescription < LCMS_Tag >::RT,
	      MZ = DimensionDescription < LCMS_Tag >::MZ
		  };
	
	    /// Default constructor
	    PickedPeakSeeder();
	
	    /// destructor
	    virtual ~PickedPeakSeeder();

	    /// Copy constructor
	    PickedPeakSeeder(const PickedPeakSeeder& rhs);
	    
	    /// Assignment operator
	    PickedPeakSeeder& operator= (const PickedPeakSeeder& rhs);
	
	    /// return next seed
	    IndexSet nextSeed() throw (NoSuccessor);
	
	    static BaseSeeder* create()
	    {
	        return new PickedPeakSeeder();
	    }
	
	    static const String getProductName()
	    {
	        return "PickedPeakSeeder";
	    }
		
		protected:
			virtual void updateMembers_();
			
	    /// Finds the neighbour of the peak denoted by @p current_mz in the previous scan
	    std::vector<double>::iterator searchInScan_(const std::vector<CoordinateType>::iterator& scan_begin, const std::vector<CoordinateType>::iterator& scan_end, CoordinateType current_mz)
	    {
	      // perform binary search to find the neighbour in rt dimension
	      // 	lower_bound finds the peak with m/z current_mz or the next larger peak if this peak does not exist.
	      std::vector<CoordinateType>::iterator insert_iter = lower_bound(scan_begin,scan_end,current_mz);
	
	      // the peak found by lower_bound does not have to be the closest one, therefore we have
	      // to check both neighbours
	      if ( insert_iter == scan_end ) // we are at the and have only one choice
	      {
	      	return --insert_iter;
	      }
	      else
	      {
	        // if the found peak is at the beginning of the spectrum,
	        // there is not much we can do.
	        if ( insert_iter == scan_begin )
	        {
	            return insert_iter;
	        }
	        else // see if the next smaller one fits better
	        {
	          double delta_mz = fabs(*insert_iter - current_mz);
	          
	          --insert_iter;
	
	          if ( fabs(*insert_iter - current_mz) < delta_mz )
	          {
	          	return insert_iter; // peak to the left is closer (in m/z dimension)
	          }
	          else
	          {
	          	return ++insert_iter;    // peak to the right is closer
	          }
	        }
	      }
	
	    } // end of searchInScan_
	
	    /// Find out which charge state belongs to this distance
	    UnsignedInt distanceToCharge_(const CoordinateType& dist);
	
	    /// Sweeps through scans and detects isotopic patterns
	    void sweep_();
	
	    /// stores the retention time of each isotopic cluster
	    TableType iso_map_;
	
	    /// Pointer to the current region
	    TableType::const_iterator curr_region_;
	
	    /// indicates whether the extender has been initialized
	    bool is_initialized_;
	
	    /// upper bound for distance between charge 1 peaks
	    CoordinateType charge1_ub_;
	    /// lower bound for distance between charge 1 peaks
	    CoordinateType charge1_lb_;
	
	    /// upper bound for distance between charge 2 peaks
	    CoordinateType charge2_ub_;
	    /// lower bound for distance between charge 2 peaks
	    CoordinateType charge2_lb_;
	
	    /// upper bound for distance between charge 3 peaks
	    CoordinateType charge3_ub_;
	    /// lower bound for distance between charge 3 peaks
	    CoordinateType charge3_lb_;
	
	    /// upper bound for distance between charge 4 peaks
	    CoordinateType charge4_ub_;
	    /// lower bound for distance between charge 4 peaks
	    CoordinateType charge4_lb_;
	
	    /// upper bound for distance between charge 5 peaks
	    CoordinateType charge5_ub_;
	    /// lower bound for distance between charge 5 peaks
	    CoordinateType charge5_lb_;
	};
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PICKEDPEAKSEEDER_H
