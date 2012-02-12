// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#ifndef OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
    @brief ThresholdMower removes all peaks below a threshold.
   	
    @htmlinclude OpenMS_ThresholdMower.parameters

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI ThresholdMower
		: public DefaultParamHandler 
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    ThresholdMower();
		/// destructor
    virtual ~ThresholdMower();
	
		/// copy constructor
		ThresholdMower(const ThresholdMower& source);
		/// assignment operator
		ThresholdMower& operator=(const ThresholdMower& source);	
		// @}

		// @name Accessors
		// @{
		///

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{			
			// sort by intensity
			spectrum.sortByIntensity();
			
			// find right position to erase
			typename SpectrumType::PeakType p;
			threshold_ = ((DoubleReal)param_.getValue("threshold"));
			p.setIntensity(threshold_);
			spectrum.erase(
										spectrum.begin(),
										lower_bound(spectrum.begin(), spectrum.end(), p, typename SpectrumType::PeakType::IntensityLess())
										);
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);
		
		//TODO reimplement DefaultParamHandler::updateMembers_()
		
	private:
		DoubleReal threshold_;
		
		// @}
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
