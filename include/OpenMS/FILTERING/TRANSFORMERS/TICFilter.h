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
#ifndef OPENMS_FILTERING_TRANSFORMERS_TICFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_TICFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

namespace OpenMS
{
  /**
  	@brief TICFilter calculates TIC
		 
		@htmlinclude OpenMS_TICFilter.parameters
  
		@ingroup SpectraFilter
  */
  class OPENMS_DLLAPI TICFilter : public FilterFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// standard constructor
    TICFilter();

    /// copy constructor
    TICFilter(const TICFilter& source);

		/// destructor
		virtual ~TICFilter();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    TICFilter& operator=(const TICFilter& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new TICFilter(); }

		///
		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::ConstIterator ConstIterator;
    	double TIC = 0;
    	//double window = (double)param_.getValue("window");
			
    	for (ConstIterator it = spectrum.begin(); it != spectrum.end();++it )
    	{
      	TIC += it->getIntensity();
   		}
    	return TIC;
		}

		///
		static const String getProductName()
		{
			return "TICFilter";
		}
		// @}

  };
}
#endif //OPENMS_FILTERING_TRANSFORMERS_TICFILTER_H

