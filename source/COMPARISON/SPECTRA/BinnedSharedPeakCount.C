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

#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>

using namespace std;

namespace OpenMS
{
	BinnedSharedPeakCount::BinnedSharedPeakCount()
	  : BinnedSpectrumCompareFunctor()
	{
			setName(BinnedSharedPeakCount::getProductName());
			defaults_.setValue("normalized", 1, "is set 1 if the similarity-measurement is normalized to the range [0,1]");
			defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, defines the distance of two PrecursorPeaks for which they are supposed to be from different peptides");
			defaultsToParam_();
	}
	
	BinnedSharedPeakCount::BinnedSharedPeakCount(const BinnedSharedPeakCount& source)
	  : BinnedSpectrumCompareFunctor(source)
	{
	}
	
	BinnedSharedPeakCount::~BinnedSharedPeakCount()
	{
	}
	
	BinnedSharedPeakCount& BinnedSharedPeakCount::operator = (const BinnedSharedPeakCount& source)
	{
		if (this != &source)
		{
	  		BinnedSpectrumCompareFunctor::operator = (source);
		}
	  	return *this;
	}

	double BinnedSharedPeakCount::operator () (const BinnedSpectrum& spec) const
	{
		return operator () (spec, spec);
	}
	
	double BinnedSharedPeakCount::operator () (const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const
	{
		if(!spec1.checkCompliance(spec2))
		{
			cout << "incompatible" << endl;
			throw BinnedSpectrumCompareFunctor::IncompatibleBinning(__FILE__, __LINE__, __PRETTY_FUNCTION__, "");
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
	  	
		double score(0), sum(0);
		UInt denominator(max(spec1.getFilledBinNumber(),spec2.getFilledBinNumber())), shared_Bins(min(spec1.getBinNumber(),spec2.getBinNumber()));
			
		// all bins at equal position that have both intensity > 0 contribute positively to score
		for (Size i = 0; i < shared_Bins; ++i)
		{			
			if(spec1.getBins()[i]>0 && spec2.getBins()[i]>0) 
			{
				sum++;
			}
		}
						
		// resulting score normalized to interval [0,1]
	    score = sum / denominator;
	
	    return score;
	
	}

}
