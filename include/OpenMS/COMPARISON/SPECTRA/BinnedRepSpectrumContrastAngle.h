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

#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDREPSPECTRUMCONTRASTANGLE_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDREPSPECTRUMCONTRASTANGLE_H

#include <OpenMS/COMPARISON/SPECTRA/BinnedRepCompareFunctor.h>

namespace OpenMS
{

  /**
  	@brief calculates the spectral contrast angle between two spectra in bin representation

		the functor does not cover the whole function, the normalization is done by the ClusterRun to allow for more flexibility.
		The SCA was defined by ???? in ????

		@ingroup SpectraComparison
  */
  class BinnedRepSpectrumContrastAngle
    : public BinnedRepCompareFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// standard constructor
    BinnedRepSpectrumContrastAngle();

    /// copy constructor
    BinnedRepSpectrumContrastAngle(const BinnedRepSpectrumContrastAngle& source);

    /// destructor
    virtual ~BinnedRepSpectrumContrastAngle();
		// @}

		// @name Operators
		// @{
		/// assignment operator
    BinnedRepSpectrumContrastAngle& operator = (const BinnedRepSpectrumContrastAngle& source);

		///
		double operator () (const BinnedRep& a, const BinnedRep& b) const;

		double operator () (const BinnedRep& a) const;
		// @}

		// @name Accessors
		// @{
		///
    static BinnedRepCompareFunctor* create(){return new BinnedRepSpectrumContrastAngle();}

		///
		static const String getProductName()
		{
			return "BinnedRepSpectrumContrastAngle";
		}
		// @}
  };

}

#endif //OPENMS_COMPARISON_SPECTRA_BINNEDREPSPECTRUMCONTRASTANGLE_H

