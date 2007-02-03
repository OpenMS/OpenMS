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
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>

#include <OpenMS/CONCEPT/Factory.h>

#include <cmath>
#include <sstream>

using namespace std;

namespace OpenMS
{
  PeakSpectrumCompareFunctor::PeakSpectrumCompareFunctor()
		: FactoryProduct(PeakSpectrumCompareFunctor::getProductName())
  {
    defaults_.setValue("filterwindow", 2.3);
		defaultsToParam_();
  }
  
  PeakSpectrumCompareFunctor::PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source)
		:	FactoryProduct(source)
  {
  }
 
	PeakSpectrumCompareFunctor::~PeakSpectrumCompareFunctor()
	{
	}
 
  PeakSpectrumCompareFunctor& PeakSpectrumCompareFunctor::operator=(const PeakSpectrumCompareFunctor& source)
  {
		if (this != &source)
		{
    	FactoryProduct::operator=(source);
		}
    return *this;
  }
 
	void PeakSpectrumCompareFunctor::registerChildren()
	{
    Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumCheapDPCorr::getProductName(), &SpectrumCheapDPCorr::create);
    Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumPrecursorComparator::getProductName(), &SpectrumPrecursorComparator::create);
		Factory<PeakSpectrumCompareFunctor>::registerProduct(ZhangSimilarityScore::getProductName(), &ZhangSimilarityScore::create);
		Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumAlignmentScore::getProductName(), &SpectrumAlignmentScore::create);
	}
	
  /**
  does a check if comparison makes sense at all<br>
  counter-examples for useful comparisons would be spectra that have a mass 
  difference bigger than the <i>filterwindow</i>, spectra with different parent 
  charge states or bin-representation spectra with different binning parameters  
  */
	/*
	double PeakSpectrumCompareFunctor::filter(const ClusterSpectrum& a, const ClusterSpectrum& b) const
  {
    double filterwindow = (double)param_.getValue("filterwindow");
    double factor = 1;
    if (a.getParentionCharge() != b.getParentionCharge())
    {
      factor = 0;
    }
    if (usebins_)
    {
      stringstream ss;
      ss << a.getBinSize() << ":" << a.getBinSpread() << " != " << b.getBinSize() << ":" << b.getBinSpread() << " ids: " << a.id() << " " << b.id();
      // its not very informative to compare spectra with different bin sizes and spread
      if (fabs(a.getBinSize() - b.getBinSize()) > 1e-8 || a.getBinSpread() != b.getBinSpread())
      {
        throw ClusterSpectrum::WrongRepresentation(__FILE__, __LINE__, __PRETTY_FUNCTION__, ss.str().c_str());
      }
    }
    if (fabs(a.getParentMass() - b.getParentMass()) > filterwindow) return 0;
    return factor;
  }
	*/

}
