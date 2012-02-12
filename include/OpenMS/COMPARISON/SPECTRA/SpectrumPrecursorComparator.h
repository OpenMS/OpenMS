// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

namespace OpenMS
{

  /**	
  	@brief SpectrumPrecursorComparator compares just the parent mass of two spectra
		 
		@htmlinclude OpenMS_SpectrumPrecursorComparator.parameters

		@ingroup SpectraComparison
  */
  class OPENMS_DLLAPI SpectrumPrecursorComparator : public PeakSpectrumCompareFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor 
    SpectrumPrecursorComparator();

    /// copy constructor 
    SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source);

    /// destructor 
    virtual ~SpectrumPrecursorComparator();
		// @}

		// @name Operators
		// @{
    /// assignment operator 
    SpectrumPrecursorComparator& operator = (const SpectrumPrecursorComparator& source);

		double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const;

		double operator () (const PeakSpectrum& a) const;
		// @}
		
		// @name Accessors
		// @{
		///
    static PeakSpectrumCompareFunctor* create() { return new SpectrumPrecursorComparator(); }

		///
		static const String getProductName()
		{
			return "SpectrumPrecursorComparator";
		}
		// @}

  };

}

#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H
