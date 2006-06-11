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
// $Id: ParentFilter.h,v 1.4 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_PARENTFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_PARENTFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

namespace OpenMS
{
  /**
  	@brief Parentfilter returns parention charge and mass<br>
  */
  class ParentFilter : public FilterFunctor
  {
  public:
    /// standard constructor
    ParentFilter();

    /// copy constructor
    ParentFilter(const ParentFilter& source);

    /// assignment operator
    ParentFilter& operator=(const ParentFilter& source);

    /// destructor
    ~ParentFilter();

    static FactoryProduct* create() { return new ParentFilter();}

    //std::vector<double> operator()(const ClusterSpectrum& cspec);

    //String info() const;

		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
	    //vector<double> result;
    	//result.push_back(0);
    	//result.push_back(0);
    	//result.push_back(0);
/*
    	try{
      	result.at(cspec.getParentionCharge()-1) = 1;
    	}
    	catch (out_of_range&)
    	{
      	cerr << "charge state not in {1,2,3} in spec " << cspec.id() << endl;
    	}
    	result.push_back(cspec.getParentMass());
    	return result;
*/
			return 0;
		}

		static const String getName()
		{
			return "ParentFilter";
		}

  private:
    //static const String info_;
  };
}
#endif //OPENMS_FILTERING_TRANSFORMERS_PARENTFILTER_H

