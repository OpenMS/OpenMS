// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_SCALER_H
#define OPENMS_FILTERING_TRANSFORMERS_SCALER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <map>

namespace OpenMS
{
  /**
  	@brief Scaler scales the peak by ranking the peaks and assigning intensity according to rank

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI Scaler
		: public DefaultParamHandler 
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    Scaler();
    /// destructor
    virtual ~Scaler();
	
		/// copy constructor
    Scaler(const Scaler& source);
    /// assignment operator
    Scaler& operator = (const Scaler& source);

		// @}

		// @name Accessors
		// @{

		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{	
			if (spectrum.empty()) return;
			
			spectrum.sortByIntensity();
			typename SpectrumType::size_type count = spectrum.size();
			++count;
			typename SpectrumType::PeakType::IntensityType last_int = 0.0; 
			typename SpectrumType::Iterator it = spectrum.end();
			do
			{
				--it;
				if (it->getIntensity()!=last_int)
				{
					--count;
				}
				last_int = it->getIntensity();
				it->setIntensity(count);
			}
			while (it!=spectrum.begin());
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

    void filterPeakMap(PeakMap& exp);
		
		//TODO reimplement DefaultParamHandler::updateMembers_() when introducing member variables

		// @}
		
  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_SCALER_H
