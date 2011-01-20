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
#ifndef OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
  	@brief NeutralLossDiffFilter returns the total intensity ob peak pairs whose m/z difference can be explained by a neutral loss
		 
		@htmlinclude OpenMS_NeutralLossDiffFilter.parameters

		@ingroup SpectraFilter
  */
  class OPENMS_DLLAPI NeutralLossDiffFilter : public FilterFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    NeutralLossDiffFilter();

    /// copy constructor
    NeutralLossDiffFilter(const NeutralLossDiffFilter& source);

		/// destructor
		virtual ~NeutralLossDiffFilter();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    NeutralLossDiffFilter& operator = (const NeutralLossDiffFilter& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new NeutralLossDiffFilter(); }

		///
		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
    	double tolerance = (double)param_.getValue("tolerance");
    	double isodiff = 0;
    	//iterate over all peaks
    	for (int i = 0; i < (int)spectrum.size(); ++i)
    	{
      	for (int j = 1; i - j >= 0; ++j)
      	{
					double pos_diff = std::fabs(spectrum[i-j].getPosition()[0] - spectrum[i].getPosition()[0]);
        	if (std::fabs(pos_diff - 18) < tolerance || std::fabs(pos_diff - 17) < tolerance) // water and ammonium
        	{
          	isodiff += spectrum[i-j].getIntensity() + spectrum[i].getIntensity();
        	}
        	else 
					{
						if (pos_diff > 18 + tolerance)
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
			return "NeutralLossDiffFilter";
		}
		// @}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H
