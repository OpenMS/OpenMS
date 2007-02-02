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
#ifndef OPENMS_FILTERING_TRANSFORMERS_KELLERQUALITY_H
#define OPENMS_FILTERING_TRANSFORMERS_KELLERQUALITY_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>

namespace OpenMS
{
	class ClusterSpectrum;
	
  /**
  	@brief KellerQuality assigns a quality measure to a spectrum.
  	
  	It is based on the linear regression formula from their 2003 Paper
  	http://www.systemsbiology.org/PDFs/Keller.Emperical%20statistical%20model.Anal%20Chem.02.pdf

		@ingroup SpectraFilter
  */	
  class KellerQuality : public FilterFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    KellerQuality();

    /// copy constructor
    KellerQuality(const KellerQuality& source);

		/// destructor
		virtual ~KellerQuality();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    KellerQuality& operator = (const KellerQuality& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new KellerQuality(); }

		///
    double operator () (const ClusterSpectrum& spec);

		///
		static const String getProductName()
		{
			return "KellerQuality";
		}
		// @}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_KELLERQUALITY_H
