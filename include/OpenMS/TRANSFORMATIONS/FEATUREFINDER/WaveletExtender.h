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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_WAVELETEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_WAVELETEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFinder.h>


#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/KernelTraits.h>

#include <OpenMS/DATASTRUCTURES/ScanIndex.h>

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
class WaveletExtender
            : public BaseExtender
{

public:

	enum DimensionId
    {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

    typedef FeaFiTraits::IntensityType IntensityType;
    typedef FeaFiTraits::CoordinateType CoordinateType;
    typedef KernelTraits::ProbabilityType ProbabilityType;

	typedef DRawDataPoint<2>::NthPositionLess<RT> RTless;
    typedef DRawDataPoint<2>::NthPositionLess<MZ> MZless;
	
	typedef IsotopeFinder<MSExperiment<DRawDataPoint<2> > >::SweepLineHash SweepLineHash;
	typedef IsotopeFinder<MSExperiment<DRawDataPoint<2> > >::DoubleList DoubleList;
	
	typedef DPeakArray<2, DRawDataPoint<2> >::iterator PeakIterator;
	
    

    /// standard constructor
    WaveletExtender();

    /// destructor
    virtual ~WaveletExtender();

    /// return next seed
    const IndexSet& extend(const UnsignedInt seed);

    /// returns an instance of this class
    static BaseExtender* create()
    {
        return new WaveletExtender();
    }

    /// returns the name of this module
    static const String getName()
    {
        return "WaveletExtender";
    }


protected:

	bool is_initialized_;
	
	IsotopeFinder<MSExperiment<DRawDataPoint<2> > >::SweepLineHash hash_;
	
	SweepLineHash::const_iterator hash_iter;
		
	ScanIndex<DPeakArray<2, DRawDataPoint<2> >  > scan_index_;
	
	CoordinateType avMZSpacing_;
	
	CoordinateType min_mass_;
	
	void copyData_(MSExperiment<DRawDataPoint<2> > & exp, DPeakArray<2, DRawDataPoint<2> >& peaks);
	
}; // end of class WaveletExtender

} // end of namespace OpenMS
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_WAVELETEXTENDER_H
