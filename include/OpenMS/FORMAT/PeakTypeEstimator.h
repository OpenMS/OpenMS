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
// $Maintainer: Marc Sturm $
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
  class PeakTypeEstimator
  {
    public:
    
    template <typename PeakConstIterator>
    SpectrumSettings::SpectrumType estimateType(const PeakConstIterator& begin, const PeakConstIterator& end) const
    {
    	double min = std::numeric_limits<double>::max();
    	double max = std::numeric_limits<double>::min();
    	double left = std::numeric_limits<double>::min();
    	double right = std::numeric_limits<double>::min();
    	
    	PeakConstIterator it = begin;
    	PeakConstIterator it2 = begin;
    	++it2;
    	for (; it2!=end; ++it,++it2)
    	{
    		min = std::min(min,it2->getPos()-it->getPos());
    		if (max < it2->getPos()-it->getPos())
    		{
    			left = it->getIntensity();
    			right = it2->getIntensity();
    			max = it2->getPos()-it->getPos();
    		}
    	}

    	//std::cout << "Min  : " << min << std::endl;
    	//std::cout << "Max  : " << max << std::endl;
    	//std::cout << "Left : " << left << std::endl;
    	//std::cout << "Right: " << right << std::endl;
    	
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
