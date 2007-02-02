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
#ifndef OPENMS_FILTERING_TRANSFORMERS_TRADSEQQUALITY_H
#define OPENMS_FILTERING_TRANSFORMERS_TRADSEQQUALITY_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

namespace OpenMS
{
	class ClusterSpectrum;
	
  /**
	  @brief TradSeqQuality returns a number > 0 if the sequest score are above a certain XCorr and above a certain deltaCN
	  
	  @param xcorr_1+ min XCorr for charge state 1
	  @param xcorr_2+ min XCorr for charge state 2
	  @param xcorr_3+ min XCorr for charge state 3
	  @param dCn_1+ min deltaCN for charge state 1
	  @param dCn_2+ min deltaCN for charge state 2
	  @param dCn_3+ min deltaCN for charge state 3

		@ingroup SpectraFilter
  */
  class TradSeqQuality : public FilterFunctor
  {
  public:
  	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    TradSeqQuality();

    /// copy constructor
    TradSeqQuality(const TradSeqQuality& source);

		/// destructor
		virtual ~TradSeqQuality();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    TradSeqQuality& operator = (const TradSeqQuality& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new TradSeqQuality(); }

		///
    double operator () (const ClusterSpectrum& spec);

		///
		static const String getProductName()
		{
			return "TradSeqQuality";
		}
		// @}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_TRADSEQQUALITY_H
