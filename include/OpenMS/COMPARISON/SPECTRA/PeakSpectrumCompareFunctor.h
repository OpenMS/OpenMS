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
#ifndef OPENMS_COMPARISON_SPECTRA_PEAKSPECTRUMCOMPAREFUNCTOR_H
#define OPENMS_COMPARISON_SPECTRA_PEAKSPECTRUMCOMPAREFUNCTOR_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{

	/**
	
		@brief Base class for compare functors of spectra; compare functors returns a similiarity value of two spectra
	
  	PeakSpectrumCompareFunctor classes return a value for a pair of PeakSpectrum objects
  	ideally the value should reflect the similarity of the pair
  	similarities of spectra should be > 0
		
		@ingroup SpectraComparison
  */
  class PeakSpectrumCompareFunctor : public FactoryProduct
  {

  public:

    /// default constructor
    PeakSpectrumCompareFunctor();

    /// copy constructor
    PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source);

    /// destructor
    virtual ~PeakSpectrumCompareFunctor();

    /// assignment operator
    PeakSpectrumCompareFunctor& operator = (const PeakSpectrumCompareFunctor& source);

    /// function call operator, calculates the similarity
    virtual double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const = 0;

		/// calculates self similarity
		virtual double operator () (const PeakSpectrum& a) const = 0;

		/// registers all derived products 
		static void registerChildren();

		/// 
		static const String getProductName()
		{
			return "PeakSpectrumCompareFunctor";
		}

  protected:

  };

}
#endif // OPENMS_COMPARISON_SPECTRA_COMPAREFUNCTOR_H
