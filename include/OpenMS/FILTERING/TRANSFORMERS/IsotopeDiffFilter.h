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
// $Id: IsotopeDiffFilter.h,v 1.4 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <string>

namespace OpenMS
{
  /**
  IsotopeDiffFilter returns total intensity of peak pairs that could result from isotope peaks<br>
  
  \param tolerance m/z tolerance
  */
  class IsotopeDiffFilter : public FilterFunctor
  {
  public:
    /** @brief standard constructor <br> */
    IsotopeDiffFilter();

    /** @brief copy constructor <br> */
    IsotopeDiffFilter(const IsotopeDiffFilter& source );

    /** @brief assignment operator <br> */
    IsotopeDiffFilter& operator=(const IsotopeDiffFilter& source );

    /** @brief destructor <br> */
    ~IsotopeDiffFilter();

    static FactoryProduct* create() { return new IsotopeDiffFilter();}

		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
   		double tolerance = (double)param_.getValue("tolerance");;
    	double isodiff = 0;
    	//iterate over all peaks
    	for (int i = 0; i < (int)spectrum.size(); ++i)
    	{
      	for (uint j = 1; i+j < spectrum.size() ; ++j)
      	{
        	if (fabs(spectrum.getContainer()[i+j].getPosition()[0] - spectrum.getContainer()[i].getPosition()[0] + 1 ) < tolerance && 
							spectrum.getContainer()[i-j].getIntensity() < spectrum.getContainer()[i].getIntensity())
        	{
          	isodiff+= spectrum.getContainer()[i].getIntensity() + spectrum.getContainer()[i+j].getIntensity();
        	}
        	else if (fabs(spectrum.getContainer()[i+j].getPosition()[0] - spectrum.getContainer()[i].getPosition()[0] ) > 1 + tolerance)
        	{
         		break;
        	}
      	}
    	}
    	//vector<double> result;
    	//result.push_back(isodiff);
    	return isodiff;
		}

		static const String getName()
		{
			return "IsotopeDiffFilter";
		}

  private:
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H
