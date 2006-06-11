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
// $Id: TICFilter.h,v 1.4 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_TICFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_TICFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

namespace OpenMS
{
  /**
  	@brief TICFilter calculates TIC and ratio of TIC in parent ion range (=parention m/z +- window )
  
  	\param window half size of parent ion range
  */
  class TICFilter : public FilterFunctor
  {
  public:
    /// standard constructor
    TICFilter();

    /// copy constructor
    TICFilter(const TICFilter& source);

    /// assignment operator
    TICFilter& operator=(const TICFilter& source);

    /// destructor
    ~TICFilter();

    static FactoryProduct* create() { return new TICFilter();}

    //std::vector<double> operator()(const ClusterSpectrum& cspec);

    //String info() const;

		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::ConstIterator ConstIterator;
			//vector<double> result;
    	double TIC = 0;
    	//double PTIC = 0;
    	double window = (double)param_.getValue("window");
			// TODO
    	//double parentpeak = cspec.getParentMass()/ cspec.getParentionCharge();
    	for (ConstIterator it = spectrum.begin(); it != spectrum.end();++it )
    	{
      	TIC += it->getIntensity();
      	//if ( fabs( it->getPosition()[0] - parentpeak )  <  window )
      	//{
        //	PTIC += it->getIntensity();
      	//}
   		}
    	//result.push_back(TIC);
    	//result.push_back(PTIC);
    	return TIC;
		}

		static const String getName()
		{
			return "TICFilter";
		}

  private:
    //static const String info_;
  };
}
#endif //OPENMS_FILTERING_TRANSFORMERS_TICFILTER_H

