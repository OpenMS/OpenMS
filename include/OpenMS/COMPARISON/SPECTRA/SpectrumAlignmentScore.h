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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <utility>

namespace OpenMS
{

  /**
	  @brief Similarity score of Zhang

		@param epsilon - defines the absolut error of the mass spectrometer; default value is 0.2 Th
		
		@ingroup SpectraComparison
  */
	
  class SpectrumAlignmentScore : public PeakSpectrumCompareFunctor
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    SpectrumAlignmentScore();

    /// copy constructor
    SpectrumAlignmentScore(const SpectrumAlignmentScore& source);

    /// destructor
    virtual ~SpectrumAlignmentScore();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    SpectrumAlignmentScore& operator = (const SpectrumAlignmentScore& source);
	
		/// 
		double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const;

		double operator () (const PeakSpectrum& spec) const;
		// @}

		// @name Accessors
		// @{
		///
    static PeakSpectrumCompareFunctor* create() { return new SpectrumAlignmentScore(); }

		///
		static const String getProductName()
		{
			return "SpectrumAlignmentScore";
		}
		// @}

  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H
