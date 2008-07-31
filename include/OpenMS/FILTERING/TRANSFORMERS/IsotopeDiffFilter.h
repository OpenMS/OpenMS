// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <string>
#include <cmath>

namespace OpenMS
{
  /**
  	@brief IsotopeDiffFilter returns total intensity of peak pairs that could result from isotope peaks
		 
		@ref IsotopeDiffFilter_Parameters are explained on a separate page.

		@ingroup SpectraFilter
  */
  class IsotopeDiffFilter : public FilterFunctor
  {
	
  public:
	
		// @name Constructors and Destrutors
		// @{
    /// default constructor
    IsotopeDiffFilter();

    /// copy constructor
    IsotopeDiffFilter(const IsotopeDiffFilter& source);

		/// destructor
		virtual ~IsotopeDiffFilter();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    IsotopeDiffFilter& operator = (const IsotopeDiffFilter& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new IsotopeDiffFilter(); }

		///
		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
   		double tolerance = (double)param_.getValue("tolerance");
    	double isodiff = 0;
			
    	//iterate over all peaks
    	for (int i = 0; i < (int)spectrum.size(); ++i)
    	{
      	for (uint j = 1; i + j < spectrum.size(); ++j)
      	{
					double pos_ij = spectrum.getContainer()[i+j].getPosition()[0];
					double pos_i = spectrum.getContainer()[i].getPosition()[0];
        	if (std::fabs(pos_ij - pos_i + 1) < tolerance)
        	{
          	isodiff += spectrum.getContainer()[i].getIntensity() + spectrum.getContainer()[i+j].getIntensity();
        	}
        	else 
					{
						if (std::fabs(spectrum.getContainer()[i+j].getPosition()[0] - spectrum.getContainer()[i].getPosition()[0]) > 1 + tolerance)
        		{
         			break;
        		}
					}
      	}
    	}
    	return isodiff;
		}

		///
		static const String getProductName()
		{
			return "IsotopeDiffFilter";
		}
		// @}

  private:
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H
