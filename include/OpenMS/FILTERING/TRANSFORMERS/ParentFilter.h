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
#ifndef OPENMS_FILTERING_TRANSFORMERS_PARENTFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_PARENTFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

namespace OpenMS
{
  /**
  	@brief Parentfilter returns parent-ion charge and mass

		// will be readded when functionality is implemented
		// ingroup SpectraFilter
  */
  class ParentFilter : public FilterFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    ///  constructor
    ParentFilter();

    /// copy constructor
    ParentFilter(const ParentFilter& source);

    /// destructor
    virtual ~ParentFilter();
		// @}
		
		// @name Operators
		// @{
		// assignment operator
		ParentFilter& operator = (const ParentFilter& source);
		// @}
		
		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new ParentFilter(); }

		///
		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
			/// @todo class needed any more? What exactly does this class? (andreas)
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
			return spectrum.getPrecursorPeak().getCharge();
		}

		///
		static const String getProductName()
		{
			return "ParentFilter";
		}
		// @}

  };
}
#endif //OPENMS_FILTERING_TRANSFORMERS_PARENTFILTER_H

