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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

namespace OpenMS
{
  /**
  	@brief ThresholdMower removes all peaks below a Threshold
   	
   	@htmlinclude OpenMS_ThresholdMower.parameters

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI ThresholdMower
    :	public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    ThresholdMower();

    /// copy constructor
    ThresholdMower(const ThresholdMower& source);

    /// destructor
    virtual ~ThresholdMower();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    ThresholdMower& operator=(const ThresholdMower& source);
		// @}

		// @name Accessors
		// @{
		///
    static PreprocessingFunctor* create() { return new ThresholdMower(); }

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{			
			// sort by intensity
			spectrum.sortByIntensity();
			
			// find right position to erase
			typename SpectrumType::PeakType p;
			p.setIntensity((double)param_.getValue("threshold"));
			spectrum.erase(
										spectrum.begin(),
										lower_bound(spectrum.begin(), spectrum.end(), p, typename SpectrumType::PeakType::IntensityLess())
										);
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);
		
		/// 
		static const String getProductName()
		{
			return "ThresholdMower";
		}
		// @}
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
