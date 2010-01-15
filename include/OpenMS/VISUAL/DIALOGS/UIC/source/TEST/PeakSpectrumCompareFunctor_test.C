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
// $Authors: Voker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Factory.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PeakSpectrumCompareFunctor, "$Id: PeakSpectrumCompareFunctor_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////

// pure interface class cannot test this

START_SECTION(PeakSpectrumCompareFunctor())
  NOT_TESTABLE
END_SECTION

START_SECTION(PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source))
  NOT_TESTABLE
END_SECTION

START_SECTION(~PeakSpectrumCompareFunctor())
	NOT_TESTABLE
END_SECTION

START_SECTION(PeakSpectrumCompareFunctor& operator = (const PeakSpectrumCompareFunctor& source))
  NOT_TESTABLE
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
  NOT_TESTABLE
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a) const)
  NOT_TESTABLE
END_SECTION

START_SECTION(static void registerChildren())
  PeakSpectrumCompareFunctor* c1 = Factory<PeakSpectrumCompareFunctor>::create("SpectrumCheapDPCorr");
  c1 = Factory<PeakSpectrumCompareFunctor>::create("SpectrumPrecursorComparator");
	TEST_EQUAL(c1->getName(), "SpectrumPrecursorComparator")
  c1 = Factory<PeakSpectrumCompareFunctor>::create("ZhangSimilarityScore");
	TEST_EQUAL(c1->getName(), "ZhangSimilarityScore")
	c1 = Factory<PeakSpectrumCompareFunctor>::create("SteinScottImproveScore");
	TEST_EQUAL(c1->getName(), "SteinScottImproveScore");
	c1 = Factory<PeakSpectrumCompareFunctor>::create("CompareFouriertransform");
	TEST_EQUAL(c1->getName(), "CompareFouriertransform")
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(PeakSpectrumCompareFunctor::getProductName(), "PeakSpectrumCompareFunctor")
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
