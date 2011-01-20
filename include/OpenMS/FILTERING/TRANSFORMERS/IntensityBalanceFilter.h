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
#ifndef OPENMS_FILTERING_TRANSFORMERS_INTENSITYBALANCEFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_INTENSITYBALANCEFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <utility>

namespace OpenMS
{
  /**
	  @brief IntensityBalanceFilter divides the m/z-range into ten regions and sums the
	         intensity in these regions. 
	  
	  The result is the intensity of the two bins with the highest intensity minus the intensity of the seven bins with lowest intensity.

		@ingroup SpectraFilter
  */
  class OPENMS_DLLAPI IntensityBalanceFilter : public FilterFunctor
  {
	
  public:
		
		// @name Constructors and Destructors
		// @{
    /// default constructor
    IntensityBalanceFilter();

    /// copy constructor
    IntensityBalanceFilter(const IntensityBalanceFilter& source);

		/// destructor
		virtual ~IntensityBalanceFilter();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    IntensityBalanceFilter& operator=(const IntensityBalanceFilter& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new IntensityBalanceFilter();}

		///
		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
			double bands = 10;
    	std::multimap<double, Size> band_intensity;
	    double parentmass = 0.0;
			if (!spectrum.getPrecursors().empty()) parentmass = spectrum.getPrecursors()[0].getMZ();
    	Size j = 0;
    	for (Size i = 0; i < bands; ++i)
    	{
      	double intensity = 0;

      	//bern 2004 says to only check between 300 and size
      	//but that seems inappropriate for small peptides (smallest is ca 450)
      	while (j < spectrum.size() && spectrum[j].getPosition()[0] < (parentmass-300)/bands*(i+1) +300)
      	{
        	intensity += spectrum[j++].getIntensity();
      	}
      	band_intensity.insert(std::make_pair(intensity,i));
    	}
    	j = 0;
    	double total_intensity = 0;
    	double twobiggest = 0;
    	double sevensmallest = 0;
    	for (std::multimap<double, Size>::reverse_iterator mmrit = band_intensity.rbegin(); mmrit != band_intensity.rend(); ++mmrit, ++j)
    	{
      	total_intensity += mmrit->first;
      	//take the two biggest
      	if (j < 2)
      	{
        	twobiggest+=mmrit->first;
      	}
      	//take the seven smallest
      	if (j > 2)
      	{
        	sevensmallest += mmrit->first;
      	}
    	}

    	return (twobiggest - sevensmallest) / total_intensity;
		}

		///
		static const String getProductName()
		{
			return "IntensityBalanceFilter";
		}
		// @}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_INTENSITYBALANCEFILTER_H
