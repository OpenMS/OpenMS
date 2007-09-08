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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Factory.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PeakSpectrumCompareFunctor, "$Id$")

/////////////////////////////////////////////////////////////

CHECK(PeakSpectrumCompareFunctor())
  // nothing to check
RESULT

CHECK(PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source))
  // nothing to check
RESULT

CHECK(~PeakSpectrumCompareFunctor())
  // nothing to check
RESULT

CHECK(PeakSpectrumCompareFunctor& operator = (const PeakSpectrumCompareFunctor& source))
  // nothing to check
RESULT

CHECK(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
  // nothing to check
RESULT

CHECK(double operator () (const PeakSpectrum& a) const)
  // nothing to check
RESULT

CHECK(static void registerChildren())
  PeakSpectrumCompareFunctor* c1 = Factory<PeakSpectrumCompareFunctor>::create("SpectrumCheapDPCorr");
  c1 = Factory<PeakSpectrumCompareFunctor>::create("SpectrumPrecursorComparator");
  c1 = Factory<PeakSpectrumCompareFunctor>::create("ZhangSimilarityScore");
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(PeakSpectrumCompareFunctor::getProductName(), "PeakSpectrumCompareFunctor")
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
