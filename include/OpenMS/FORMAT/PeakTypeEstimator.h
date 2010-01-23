// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEAKTYPEESTIMATOR_H
#define OPENMS_FORMAT_PEAKTYPEESTIMATOR_H

#include <OpenMS/METADATA/SpectrumSettings.h>

#include <cmath>
#include <limits>

namespace OpenMS
{
	/**
		@brief Estimates if the data of a spectrum is raw data or peak data
  
  	@ingroup Format
	*/
  class OPENMS_DLLAPI PeakTypeEstimator
  {
    public:
    	/**
    		@brief Estimates the peak type of the peaks in the iterator range
    		
    		@note if there are fewer than 5 peaks in the iterator range SpectrumSettings::UNKOWN is returned
    	*/
	    template <typename PeakConstIterator>
	    SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end) const
	    {
	    	//abort if there are less than 5 peak in the iterator range
	    	if (end - begin < 5)
	    	{
	    		return SpectrumSettings::UNKNOWN;
	    	}
	    	
	    	DoubleReal min = std::numeric_limits<DoubleReal>::max();
	    	DoubleReal max = std::numeric_limits<DoubleReal>::min();
	    	DoubleReal left = std::numeric_limits<DoubleReal>::min();
	    	DoubleReal right = std::numeric_limits<DoubleReal>::min();
	    	
	    	PeakConstIterator it = begin;
	    	PeakConstIterator it2 = begin;
	    	++it2;
	    	for (; it2!=end; ++it,++it2)
	    	{
	    		min = std::min(min,it2->getMZ()-it->getMZ());
	    		if (max < it2->getMZ()-it->getMZ())
	    		{
	    			left = it->getIntensity();
	    			right = it2->getIntensity();
	    			max = it2->getMZ()-it->getMZ();
	    		}
	    	}
	
	    	//std::cout << "Min  : " << min << "\n";
	    	//std::cout << "Max  : " << max << "\n";
	    	//std::cout << "Left : " << left << "\n";
	    	//std::cout << "Right: " << right << "\n";
	    	
	    	//raw data with zeros
	    	if ((max-min)<0.5)
	    	{
	    		return SpectrumSettings::RAWDATA;
	    	}
	    	//raw data without zeros
	    	else if (left < 2.0 && right < 2.0)
	    	{
	    		return SpectrumSettings::RAWDATA;
	    	}
	    	
	    	return SpectrumSettings::PEAKS;
	    }
			
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEAKTYPEESTIMATOR_H
