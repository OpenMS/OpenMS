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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmWavelet.h>

///////////////////////////

START_TEST(FeatureFinderAlgorithmWavelet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmWavelet<Peak1D,Feature> FFAW;

FFAW* ptr;
START_SECTION(FeatureFinderAlgorithmWavelet())
	ptr = new FFAW;
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION((virtual ~FeatureFinderAlgorithmWavelet()))
	delete ptr;
END_SECTION

//TODO: this should work (Clemens, Rene)
//START_SECTION([EXTRA] FeatureFinderAlgorithmWavelet() - with RichPeak1D)
//	FeatureFinderAlgorithmWavelet<RichPeak1D,Feature> ffa;
//END_SECTION

START_SECTION(virtual void run())
	// dummy subtest
	TEST_EQUAL(1,1)
END_SECTION

START_SECTION((virtual Param getDefaultParameters() const))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
END_SECTION

START_SECTION((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
	FeatureFinderAlgorithm<Peak1D,Feature>* ptr2 = FFAW::create();
	TEST_NOT_EQUAL(ptr2,0)
	delete ptr2;
END_SECTION

START_SECTION(static const String getProductName())
    TEST_EQUAL(FFAW::getProductName(),"wavelet")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
