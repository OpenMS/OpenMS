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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PPEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PPEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>

#include <OpenMS/DATASTRUCTURES/RunningAveragePosition.h>

#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/KernelTraits.h>

#include <OpenMS/MATH/MISC/LinearInterpolation.h>

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
 
#include <fstream>
#include<sstream>

namespace OpenMS
{

/**
  @brief Implements the extension phase of the FeatureFinder.
  
  @ingroup FeatureFinder
	
 */
class PPExtender
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
	
	typedef DRawDataPoint<1> RawDataPointType;
	typedef DPeakArrayNonPolymorphic<1, RawDataPointType > RawDataArrayType;
	 typedef RawDataArrayType::iterator RawDataPointIterator;
	
    /// standard constructor
    PPExtender();

    /// destructor
    virtual ~PPExtender();

    /// return next seed
    const IndexSet& extend(const UnsignedInt seed);

    /// returns an instance of this class
    static BaseExtender* create()
    {
        return new PPExtender();
    }

    /// returns the name of this module
    static const String getName()
    {
        return "PPExtender";
    }

    /// stores information about an isotopic cluser (i.e. potential peptide charge variant)
    struct IsotopeCluster
    {
        IsotopeCluster()
          : charge_(0), peaks_(), scans_()
        {}
     
        // predicted charge state of this peptide
        UnsignedInt charge_;
        // peaks in this cluster
        std::vector< UnsignedInt > peaks_;
        // the scans of this cluster
        std::vector<CoordinateType> scans_;
    };
		
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
	
	
	void getMaxPositions_( RawDataPointIterator first, RawDataPointIterator last, const ContinuousWaveletTransform& wt, std::vector<int>& localmax, unsigned int curr_peak);
 
  	/// Sums a scan
  	void sumUp_(RawDataArrayType& scan, unsigned int current_index);
	
	void AlignAndSum_(RawDataArrayType& scan, RawDataArrayType& neighbour);
 	
	/// Test if the distance between two peaks is equal to 1  / z (where z=1,2,....)
	UnsignedInt testDistance2NextPeak_(CoordinateType dist2nextpeak);

    /// Sweeps through scans and detects isotopic patterns
    void sweep_();

    /// stores the retention time of each isotopic cluster
    std::map<CoordinateType,IsotopeCluster> iso_map_;

    /// Pointer to the current region
    std::map<CoordinateType,IsotopeCluster>::const_iterator curr_region_;

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
	
	/// Computes the wavelet transform for a given scan
	ContinuousWaveletTransformNumIntegration cwt_;
	
	/// wavelet dilation
	float cwt_scale_;
				
	/// Minimum ion count
	double noise_level_signal_;

    /// The min. intensity in the cwt 
    float noise_level_cwt_;


	
}; // end of class PPExtender

} // end of namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PPEXTENDER_H
