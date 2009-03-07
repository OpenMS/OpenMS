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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <cmath>

using namespace std;

namespace OpenMS
{
  SpectrumPrecursorComparator::SpectrumPrecursorComparator()
    : PeakSpectrumCompareFunctor()
  {
		setName(SpectrumPrecursorComparator::getProductName());
    defaults_.setValue("window", 2, "Allowed deviation between precursor peaks.");
		defaultsToParam_();
  }

  SpectrumPrecursorComparator::SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source)
    : PeakSpectrumCompareFunctor(source)
  {
  }

  SpectrumPrecursorComparator::~SpectrumPrecursorComparator()
  {
  }

  SpectrumPrecursorComparator& SpectrumPrecursorComparator::operator=(const SpectrumPrecursorComparator& source)
  {
		if (this != &source)
		{
    	PeakSpectrumCompareFunctor::operator=(source);
		}
    return *this;
  }

	double SpectrumPrecursorComparator::operator () (const PeakSpectrum& spec) const
	{
		return operator () (spec, spec);
	}

  double SpectrumPrecursorComparator::operator()(const PeakSpectrum& x, const PeakSpectrum& y)const
  {
    double window = (double)param_.getValue("window");
    
    DoubleReal mz1 = 0.0;
    if (!x.getPrecursors().empty()) mz1 = x.getPrecursors()[0].getMZ();
    DoubleReal mz2 = 0.0;
    if (!y.getPrecursors().empty()) mz2 = y.getPrecursors()[0].getMZ();
    
    if (fabs (mz1 - mz2) > window) 
		{
			return 0;
		}

    return window - fabs(mz1 - mz2);
  }

}

