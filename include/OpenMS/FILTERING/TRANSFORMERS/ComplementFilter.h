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
// $Id: ComplementFilter.h,v 1.4 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>

namespace OpenMS
{
  /**
  total intensity of peak pairs that could result from complementing fragments of charge state 1<br>
  
  \param tolerance mass tolerance <br>
  */
  class ComplementFilter : public FilterFunctor
  {
  public:
    /** @brief standard constructor <br> */
    ComplementFilter();

    /** @brief copy constructor <br> */
    ComplementFilter(const ComplementFilter& source );

    /** @brief assignment operator <br> */
    ComplementFilter& operator=(const ComplementFilter& source );

    /** @brief destructor <br> */
    ~ComplementFilter();

    static FactoryProduct* create() { return new ComplementFilter();}

    //std::vector<double> operator()(const ClusterSpectrum& spec);

		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
	    double tolerance = (double)param_.getValue("tolerance");
	    double parentmass = spectrum.getPrecursorPeak().getPosition()[0];
	    double result(0);
			
	    uint j = spectrum.size() - 1;
	    for (uint i = 0; i < spectrum.size() && i <= j; ++i)
	    {
	      double sum = spectrum.getContainer()[i].getPosition()[0] + spectrum.getContainer()[j].getPosition()[0];
	      if (fabs(sum - parentmass) < tolerance )
	      {
	        result += spectrum.getContainer()[i].getIntensity() + spectrum.getContainer()[j].getIntensity();
	      }
				else if (sum < parentmass)
	      {
	        ++i;
	      }
	      else if (sum > parentmass)
	      {
	        --j;
	      }
	      else
	      {
	        // @todo another Exception?!?  Which one? ????
	        //throw CanNotHappen(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	      }
	    }
	    //vector<double> res;
	    //res.push_back(result);
	    return result;
		}

		static const String getName()
		{
			return "ComplementFilter";
		}

  private:

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTFILTER_H
