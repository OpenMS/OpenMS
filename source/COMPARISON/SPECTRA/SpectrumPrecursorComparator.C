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
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  SpectrumPrecursorComparator::SpectrumPrecursorComparator()
    : PeakSpectrumCompareFunctor()
  {
		setName(SpectrumPrecursorComparator::getProductName());
    defaults_.setValue("window", 2);
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
    //double filterfactor = filter(csa,csb);
    double score = 0;
    double window = (double)param_.getValue("window");
    //const MSSpectrum< DPeak<1> >& x = csa.getSpec();
    //const MSSpectrum< DPeak<1> >& y = csb.getSpec();
    
    if (fabs (x.getPrecursorPeak().getPosition()[0] - y.getPrecursorPeak().getPosition()[0]) > window) 
		{
			return 0;
		}
    score = window - fabs(x.getPrecursorPeak().getPosition()[0] - y.getPrecursorPeak().getPosition()[0]);
    return score;
  }

}

