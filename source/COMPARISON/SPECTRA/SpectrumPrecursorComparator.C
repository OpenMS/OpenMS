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
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  SpectrumPrecursorComparator::SpectrumPrecursorComparator()
    : CompareFunctor()
  {
		name_ = SpectrumPrecursorComparator::getName();
    defaults_.setValue("window", 2);
		param_ = defaults_;
    usebins_ = false;
  }

  SpectrumPrecursorComparator::SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source)
    : CompareFunctor(source)
  {
  }

  SpectrumPrecursorComparator::~SpectrumPrecursorComparator()
  {
  }

  SpectrumPrecursorComparator& SpectrumPrecursorComparator::operator=(const SpectrumPrecursorComparator& source)
  {
    CompareFunctor::operator=(source);
    return *this;
  }


  double SpectrumPrecursorComparator::operator()(const ClusterSpectrum& csa, const ClusterSpectrum& csb)const
  {
    double filterfactor = filter(csa,csb);
    double score = 0;
    double window = (double)param_.getValue("window");
    const MSSpectrum< DPeak<1> >& x = csa.getSpec();
    const MSSpectrum< DPeak<1> >& y = csb.getSpec();
    
    if ( fabs (x.getPrecursorPeak().getPosition()[0] - y.getPrecursorPeak().getPosition()[0]) > window ) return 0;
    score = window - fabs (x.getPrecursorPeak().getPosition()[0] - y.getPrecursorPeak().getPosition()[0]);
    return filterfactor * score;
  }

}

