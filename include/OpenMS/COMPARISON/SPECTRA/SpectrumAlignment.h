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
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <utility>

namespace OpenMS
{

  /**
	  @brief Aligns the peaks of two spectra

		@param epsilon - defines the absolut error of the mass spectrometer; default value is 0.3
		
		@ingroup SpectraComparison
  */
	
  class SpectrumAlignment : public FactoryProduct
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    SpectrumAlignment();

    /// copy constructor
    SpectrumAlignment(const SpectrumAlignment& source);

    /// destructor
    virtual ~SpectrumAlignment();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    SpectrumAlignment& operator = (const SpectrumAlignment& source);
		// @}

		// @name Accessors
		// @{
		///
		
		///
		static const String getProductName()
		{
			return "SpectrumAlignment";
		}

		//template <typename SpectrumType> void getSpectrumAlignment(std::vector<std::pair<Size, Size> >& alignment, const SpectrumType& s1, const SpectrumType& s2) const;

		// TODO code from this method into template method above and this method should simply call the method above
		void getSpectrumAlignment(std::vector<std::pair<Size, Size> >& alignment, const PeakSpectrum& s1, const PeakSpectrum& s2) const;
		// @}

  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H
