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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <string>
#include <cmath>

namespace OpenMS
{
  /**
  	@brief GoodDiffFilter counts the number ob peak pairs whose m/z difference can be explained by a amino acid loss
		 
		@ref GoodDiffFilter_Parameters are explained on a separate page.

		@ingroup SpectraFilter
  */
  class GoodDiffFilter : public FilterFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    GoodDiffFilter();

    /// copy constructor
    GoodDiffFilter(const GoodDiffFilter& source);

		/// destructor
		virtual ~GoodDiffFilter();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    GoodDiffFilter& operator=(const GoodDiffFilter& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new GoodDiffFilter();}

		///
		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
	    double tolerance = (double)param_.getValue("tolerance");
  	  double gooddiff = 0;
    	//iterate over all peaks
    	double totaldiff = 0;
    	for (uint i = 0; i < spectrum.size(); ++i)
    	{
      	//look for each peakdifference that is in range of aa residuemasses (56/187), if it could be a aa (aamass)
      	for (uint j = i; i+j < spectrum.size(); ++j)
      	{
        	double diff =  spectrum.getContainer()[i+j].getPosition()[0] - spectrum.getContainer()[i].getPosition()[0];
        	if (diff < 56)
        	{
						continue;
					}
					
					if (diff > 187)
        	{
          	j = spectrum.size();
        	}
        	else
        	{
          	totaldiff += spectrum.getContainer()[i+j].getIntensity() + spectrum.getContainer()[i].getIntensity();
          	std::map<double, char>::const_iterator aait = aamass_.lower_bound(diff);
          	//look for aamasses that fit diff
          	if (fabs(aait->first - diff ) <= tolerance)
          	{
           		gooddiff += spectrum.getContainer()[i+j].getIntensity()  + spectrum.getContainer()[i].getIntensity();
          	}
          	else
          	{
           		++aait;
           		if ((aait) != aamass_.end() && fabs ((aait)->first - diff) <= tolerance)
           		{
             		gooddiff += spectrum.getContainer()[i+j].getIntensity() + spectrum.getContainer()[i].getIntensity();
           		}
						}
          }
      	}
    	}

    	return gooddiff/totaldiff;
		}

		///
		static const String getProductName()
		{
			return "GoodDiffFilter";
		}
		// @}


		private: 
		
    	/// list of unique amino acid masses
    	std::map<double, char> aamass_;
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H
