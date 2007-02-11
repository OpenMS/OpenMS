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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>

namespace OpenMS
{
	/** 
		@brief Seeding module which selects single peaks based on their s/n ratio.
		
		Groups of peaks a clustered within a certain distance and traced over consecutive scans.
		
		This class we developed for a third-party project and we should discuss whether we want
		it to become part of OpenMS.
		
		@ingroup FeatureFinder
	*/ 
  class PeakSeeder 
    : public BaseSeeder
  {
  	public:	
			typedef FeaFiTraits::IntensityType IntensityType;
		  typedef FeaFiTraits::CoordinateType CoordinateType;
		  typedef KernelTraits::ProbabilityType ProbabilityType;	
	
	    enum DimensionId
	    {
	        RT = DimensionDescription < LCMS_Tag >::RT,
	        MZ = DimensionDescription < LCMS_Tag >::MZ
	    };
	
			typedef FeaFiTraits::MapType MapType;
			typedef MapType::SpectrumType SpectrumType;
			typedef MapType::PeakType PeakType;
			typedef std::multimap<CoordinateType,IsotopeCluster> TableType;
	
		  /// Default constructor
	    PeakSeeder();
	
	    /// destructor 
	    virtual ~PeakSeeder();

	    /// Copy constructor
	    PeakSeeder(const PeakSeeder& rhs);
	    
	    /// Assignment operator
	    PeakSeeder& operator= (const PeakSeeder& rhs);
	
	    /// return next seed 
	    IndexSet nextSeed() throw (NoSuccessor);
	
	    static BaseSeeder* create()
	    {
	      return new PeakSeeder();
	    }
	
	    static const String getProductName()
	    {
	      return "PeakSeeder";
	    }
	
	  protected:
	 
	    /// Finds the neighbour of the peak denoted by @p current_mz in the previous scan
	    std::vector<double>::const_iterator searchInScan_(const std::vector<CoordinateType>::const_iterator& scan_begin, const std::vector<CoordinateType>::const_iterator& scan_end, CoordinateType current_mz)
	    {
	      // perform binary search to find the neighbour in rt dimension
	      // 	lower_bound finds the peak with m/z current_mz or the next larger peak if this peak does not exist.
	      std::vector<CoordinateType>::const_iterator insert_iter = lower_bound(scan_begin,scan_end,current_mz);
	
	      // the peak found by lower_bound does not have to be the closest one, therefore we have
	      // to check both neighbours
	      if ( insert_iter == scan_end ) // we are at the and have only one choice
	      {
	      	--insert_iter;
	      }
        // if the found peak is at the beginning of the spectrum,
        // there is not much we can do.
        else if ( insert_iter != scan_begin )
        {
          if ( *insert_iter - current_mz < current_mz - *(--insert_iter) )
          {
          	++insert_iter;    // peak to the right is closer
          }
	      }

				return insert_iter;
	    }

			/// Finds local maxima in the cwt
			void filterAndComputeLocalMax_(const SpectrumType & vec, 
														 						 							std::vector<int>& localmax
																											#ifdef DEBUG_FEATUREFINDER
														 													,const UnsignedInt& currscan_index
																											#endif
														 													);
																											
			/// Retrieves the iterator for a peak cluster at mz @p  curr_mz.
			TableType::iterator retrieveHashIter_(const CoordinateType& curr_mz, 
																													 CoordinateType& mz_in_hash, 
																													 const std::vector<CoordinateType>& iso_last_scan,
																													 const UnsignedInt& currscan_index );
	
		  /// Sweeps through scans and detects isotopic patterns
		  void sweep_();
		
		  /// stores the retention time of each isotopic cluster
		  TableType iso_map_;
		
		  /// Pointer to the current region
		  TableType::const_iterator curr_region_;
		
		  /// indicates whether the extender has been initialized
		  bool is_initialized_;
			
			/// 
			
	 
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKSEEDER_H
