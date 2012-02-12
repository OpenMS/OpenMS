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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSHAREDPEAKCOUNT_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSHAREDPEAKCOUNT_H

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>

#include <cfloat>
#include <cmath>

namespace OpenMS
{

  /**
	  @brief Compare functor scoring the shared peaks for similarity measurement

		The details of the score can be found in:
		K. Wan, I. Vidavsky, and M. Gross. Comparing similar spectra: from
		similarity index to spectral contrast angle. Journal of the American Society
		for Mass Spectrometry, 13(1):85{88, January 2002.

		@htmlinclude OpenMS_BinnedSharedPeakCount.parameters

		@see BinnedSpectrumCompareFunctor @see BinnedSpectrum

		@ingroup SpectraComparison
  */

  class OPENMS_DLLAPI BinnedSharedPeakCount : public BinnedSpectrumCompareFunctor
  {
  public:

    /// default constructor
    BinnedSharedPeakCount();

    /// copy constructor
    BinnedSharedPeakCount(const BinnedSharedPeakCount& source);

    /// destructor
    virtual ~BinnedSharedPeakCount();

    /// assignment operator
    BinnedSharedPeakCount& operator = (const BinnedSharedPeakCount& source);

    /** function call operator, calculates the similarity of the given arguments

				@param spec1 First spectrum given as a binned representation
				@param spec2 Second spectrum given as a binned representation
				@throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same
		*/
		double operator () (const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const;

	/// function call operator, calculates self similarity
	double operator () (const BinnedSpectrum& spec) const;

	///
    static BinnedSpectrumCompareFunctor* create() { return new BinnedSharedPeakCount(); }

	/// get the identifier for this DefaultParamHandler
	static const String getProductName()
	{
		return "BinnedSharedPeakCount";
	}

  };

}
#endif //OPENMS_COMPARISON_SPECTRA_BINNEDSHAREDPEAKCOUNT_H
