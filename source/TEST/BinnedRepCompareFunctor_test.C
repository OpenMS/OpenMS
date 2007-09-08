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

#include <OpenMS/COMPARISON/SPECTRA/BinnedRepCompareFunctor.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(BinnedRepCompareFunctor, "$Id$")

/////////////////////////////////////////////////////////////

CHECK(BinnedRepCompareFunctor())
	// nothing to check
RESULT

CHECK(BinnedRepCompareFunctor(const BinnedRepCompareFunctor& source))
	// nothing to check
RESULT

CHECK(~BinnedRepCompareFunctor())
	// nothing to check
RESULT

CHECK(BinnedRepCompareFunctor& operator = (const BinnedRepCompareFunctor& source))
	// nothing to check
RESULT

CHECK(double operator () (const BinnedRep& s1, const BinnedRep& s2) const)
	// nothing to check
RESULT

CHECK(double operator () (const BinnedRep& a) const)
	// nothing to check
RESULT

CHECK(static void registerChildren())
	BinnedRepCompareFunctor* c1 = Factory<BinnedRepCompareFunctor>::create("BinnedRepSpectrumContrastAngle");
	c1 = Factory<BinnedRepCompareFunctor>::create("BinnedRepSharedPeakCount");
	c1 = Factory<BinnedRepCompareFunctor>::create("BinnedRepSumAgreeingIntensities");
	c1 = Factory<BinnedRepCompareFunctor>::create("BinnedRepMutualInformation");
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(BinnedRepCompareFunctor::getProductName(), "BinnedRepCompareFunctor")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
