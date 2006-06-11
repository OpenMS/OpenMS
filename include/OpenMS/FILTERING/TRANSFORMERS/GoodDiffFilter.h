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
// $Id: GoodDiffFilter.h,v 1.4 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <string>

namespace OpenMS
{
  /**
  GoodDiffFilter counts the number ob peak pairs whose m/z difference can be explained by a amino acid loss <br>
  
  \param tolerance m/z tolerance
  */
  class GoodDiffFilter : public FilterFunctor
  {
  public:
    /** @brief standard constructor <br> */
    GoodDiffFilter();

    /** @brief copy constructor <br> */
    GoodDiffFilter(const GoodDiffFilter& source);

    /** @brief assignment operator <br> */
    GoodDiffFilter& operator=(const GoodDiffFilter& source);

    /** @brief destructor <br> */
    ~GoodDiffFilter();

    static FactoryProduct* create() { return new GoodDiffFilter();}

    //std::vector<double> operator()(const ClusterSpectrum& spec);

    //String info() const;
	
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
        	{}
        	else if (diff > 187)
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
    	//vector<double> result;
    	//result.push_back(gooddiff/totaldiff);
    	return gooddiff/totaldiff;
		}


		static const String getName()
		{
			return "GoodDiffFilter";
		}

  private:
    //static const String info_;

    /**
    list of unique amino acid masses<br>
    */
    std::map<double, char> aamass_;
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H
