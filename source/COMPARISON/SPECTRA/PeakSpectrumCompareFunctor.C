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
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SteinScottImproveScore.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakAlignment.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <cmath>
#include <sstream>

using namespace std;

namespace OpenMS
{
  PeakSpectrumCompareFunctor::PeakSpectrumCompareFunctor()
		: DefaultParamHandler(PeakSpectrumCompareFunctor::getProductName())
  {
  }
  
  PeakSpectrumCompareFunctor::PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source)
		:	DefaultParamHandler(source)
  {
  }
 
	PeakSpectrumCompareFunctor::~PeakSpectrumCompareFunctor()
	{
	}
 
  PeakSpectrumCompareFunctor& PeakSpectrumCompareFunctor::operator=(const PeakSpectrumCompareFunctor& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator=(source);
		}
    return *this;
  }
 
	void PeakSpectrumCompareFunctor::registerChildren()
	{
		Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumCheapDPCorr::getProductName(), &SpectrumCheapDPCorr::create);
    Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumPrecursorComparator::getProductName(), &SpectrumPrecursorComparator::create);
		Factory<PeakSpectrumCompareFunctor>::registerProduct(ZhangSimilarityScore::getProductName(), &ZhangSimilarityScore::create);
		Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumAlignmentScore::getProductName(), &SpectrumAlignmentScore::create);
		Factory<PeakSpectrumCompareFunctor>::registerProduct(SteinScottImproveScore::getProductName(), &SteinScottImproveScore::create);
		Factory<PeakSpectrumCompareFunctor>::registerProduct(CompareFouriertransform::getProductName(), &CompareFouriertransform::create);
		Factory<PeakSpectrumCompareFunctor>::registerProduct(PeakAlignment::getProductName(), &PeakAlignment::create);
	}

}
