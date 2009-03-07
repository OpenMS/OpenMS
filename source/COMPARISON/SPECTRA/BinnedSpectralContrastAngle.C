// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>

using namespace std;

namespace OpenMS
{
	BinnedSpectralContrastAngle::BinnedSpectralContrastAngle()
	  : BinnedSpectrumCompareFunctor()
	{
			setName(BinnedSpectralContrastAngle::getProductName());
			defaults_.setValue("normalized", 1, "is set 1 if the similarity-measurement is normalized to the range [0,1]");
			defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, defines the distance of two PrecursorPeaks for which they are supposed to be from different peptides");
			defaultsToParam_();
	}
	
	BinnedSpectralContrastAngle::BinnedSpectralContrastAngle(const BinnedSpectralContrastAngle& source)
	  : BinnedSpectrumCompareFunctor(source)
	{
	}
	
	BinnedSpectralContrastAngle::~BinnedSpectralContrastAngle()
	{
	}
	
	BinnedSpectralContrastAngle& BinnedSpectralContrastAngle::operator = (const BinnedSpectralContrastAngle& source)
	{
		if (this != &source)
		{
	  		BinnedSpectrumCompareFunctor::operator = (source);
		}
	  	return *this;
	}

	double BinnedSpectralContrastAngle::operator () (const BinnedSpectrum& spec) const
	{
		return operator () (spec, spec);
	}
	
	double BinnedSpectralContrastAngle::operator () (const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const
	{
		if(!spec1.checkCompliance(spec2))
		{
			throw IncompatibleBinning(__FILE__, __LINE__, __PRETTY_FUNCTION__, "");
		}
		
		// shortcut similarity calculation by comparing PrecursorPeaks (PrecursorPeaks more than delta away from each other are supposed to be from another peptide)
		DoubleReal pre_mz1 = 0.0;
		if (!spec1.getPrecursors().empty()) pre_mz1 = spec1.getPrecursors()[0].getMZ();
		DoubleReal pre_mz2 = 0.0;
		if (!spec1.getPrecursors().empty()) pre_mz2 = spec2.getPrecursors()[0].getMZ();
		if(fabs(pre_mz1-pre_mz2)>(double)param_.getValue("precursor_mass_tolerance"))
		{
			return 0;
		}
	  			
		double score(0), numerator(0), sharedBins(min(spec1.getBinNumber(),spec2.getBinNumber())), sum1(0), sum2(0);
				
		// all bins at equal position that have both intensity > 0 contribute positively to score
		for (Size i = 0; i < sharedBins; ++i)
		{
			sum1 += spec1.getBins()[i]*spec1.getBins()[i];
			sum2 += spec2.getBins()[i]*spec2.getBins()[i];
			numerator += (spec1.getBins()[i] * spec2.getBins()[i]);
		}
		
		// resulting score standardized to interval [0,1]
		score = numerator / (sqrt(sum1 * sum2));
		
		return score;
	
	}

}
