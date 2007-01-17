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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MARRWAVELETSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MARRWAVELETSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>

#include <OpenMS/DATASTRUCTURES/IndexSet.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/KernelTraits.h>

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
 
#include <fstream>
#include<sstream>

namespace OpenMS
{

	/** @brief Seeding class 
	  
	 		Uses the continuous wavelet transform (Marr mother wavelet) to detect isotopic pattern
			in each scan. Patterns that occur in several consecutive scans are declared as seeds
			for the extension phase.
				
			@ingroup FeatureFinder
		
	*/ 
  class MarrWaveletSeeder 
    : public BaseSeeder
  {
	
	typedef FeaFiTraits::IntensityType IntensityType;
  typedef FeaFiTraits::CoordinateType CoordinateType;
  typedef KernelTraits::ProbabilityType ProbabilityType;
	
	/// stores information about an isotopic cluser (i.e. potential peptide charge variant)
    struct IsotopeCluster
    {
        IsotopeCluster()
          : charge_(0), peaks_(), scans_()
        {}
     
        // predicted charge state of this peptide
        UnsignedInt charge_;
        // peaks in this cluster
        IndexSet peaks_;
        // the scans of this cluster
        std::vector<CoordinateType> scans_;
    };
	
	

    enum DimensionId
    {
        RT = DimensionDescription < LCMS_Tag >::RT,
        MZ = DimensionDescription < LCMS_Tag >::MZ
    };

	typedef FeaFiTraits::MapType MapType;
	typedef MapType::SpectrumType SpectrumType;
	typedef SpectrumType::ContainerType ContainerType;
	typedef	 MapType::PeakType PeakType;
	typedef SpectrumType::iterator RawDataPointIterator;
	typedef std::multimap<CoordinateType,IsotopeCluster> TableType;
	typedef TableType::value_type TableEntry;
	
  public:
	   /// standard constructor
    MarrWaveletSeeder();

    /// destructor 
    virtual ~MarrWaveletSeeder();

    /// return next seed 
    IndexSet nextSeed() throw (NoSuccessor);

    static BaseSeeder* create()
    {
      return new MarrWaveletSeeder();
    }

    static const String getName()
    {
      return "MarrWaveletSeeder";
    }
		
		

  protected:
  	/// Finds the neighbour of the peak denoted by @p current_mz in the previous scan
    std::vector<double>::iterator searchInScan_(std::vector<CoordinateType>::iterator scan_begin,
            																												std::vector<CoordinateType>::iterator scan_end ,
            																												double current_mz)
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
	
	/// Finds local maxima in the cwt
	void getMaxPositions_( RawDataPointIterator first, RawDataPointIterator last, const ContinuousWaveletTransform& wt, std::vector<int>& localmax,CoordinateType curr_peak);
 
  /// Sums the intensities in adjacent scans
  void sumUp_(ContainerType& scan, UnsignedInt& current_scan_index );
	
	/// Aligns to scans
	void AlignAndSum_(ContainerType& scan, ContainerType& neighbour);
 	
	/// Test if the distance between two peaks is equal to 1  / z (where z=1,2,....)
	UnsignedInt testDistance2NextPeak_(CoordinateType dist2nextpeak);

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
	/// lower bound for distance between charge 4 peaks
	CoordinateType charge5_lb_;
	
	/// Computes the wavelet transform for a given scan
	ContinuousWaveletTransformNumIntegration cwt_;
				
	/// Minimum ion count
	IntensityType noise_level_signal_;

   /// The min. intensity in the cwt 
   IntensityType noise_level_cwt_;

	 /// S/N threshold for single peaks in the cwt
	 IntensityType high_peak_intensity_factor_;
	 
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MARRWAVELETSEEDER_H
