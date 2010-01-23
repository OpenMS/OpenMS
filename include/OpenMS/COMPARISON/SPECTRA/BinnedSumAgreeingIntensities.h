// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSUMAGREEINGINTENSITIES_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSUMAGREEINGINTENSITIES_H

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>

#include <cmath>
#include <cfloat>

namespace OpenMS
{

  /**
	  @brief Compare functor scoring the sum of agreeing intensities for similarity measurement

		Transformation and other factors of the peptide mass spectrometry pairwise peak-list comparison process
		Witold E Wolski , Maciej Lalowski* , Peter Martus* , Ralf Herwig* , Patrick Giavalisco , Johan Gobom , Albert Sickmann , Hans Lehrach and Knut Reinert*
		BMC Bioinformatics 2005, 6:285     doi:10.1186/1471-2105-6-285

		@htmlinclude OpenMS_BinnedSumAgreeingIntensities.parameters

		@see BinnedSpectrumCompareFunctor @see BinnedSpectrum

		@ingroup SpectraComparison
  */

  class OPENMS_DLLAPI BinnedSumAgreeingIntensities : public BinnedSpectrumCompareFunctor
  {
  public:

    /// default constructor
    BinnedSumAgreeingIntensities();

    /// copy constructor
    BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities& source);

    /// destructor
    virtual ~BinnedSumAgreeingIntensities();

    /// assignment operator
    BinnedSumAgreeingIntensities& operator = (const BinnedSumAgreeingIntensities& source);

    /** function call operator, calculates the similarity of the given arguments

				@param spec1 First spectrum given as a binned representation
				@param spec2 Second spectrum given as a binned representation
				@throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same
		*/
		double operator () (const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const;

		/// function call operator, calculates self similarity
		double operator () (const BinnedSpectrum& spec) const;

	///
    static BinnedSpectrumCompareFunctor* create() { return new BinnedSumAgreeingIntensities(); }

	/// get the identifier for this DefaultParamHandler
	static const String getProductName()
	{
		return "BinnedSumAgreeingIntensities";
	}

  };

}
#endif //OPENMS_COMPARISON_SPECTRA_BINNEDSUMAGREEINGINTENSITIES_H
