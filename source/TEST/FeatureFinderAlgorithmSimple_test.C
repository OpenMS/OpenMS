// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Marcel Grunert $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSimple.h>

///////////////////////////

START_TEST(FeatureFinderAlgorithmSimple, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmSimple<Peak1D,Feature> FFASS;

FFASS* ptr;
CHECK(FeatureFinderAlgorithmSimple())
	ptr = new FFASS;
	TEST_NOT_EQUAL(ptr,0)
	delete ptr;
RESULT

CHECK(virtual void run())
	// dummy subtest
	TEST_EQUAL(1,1)
RESULT

CHECK((virtual Param getDefaultParameters() const))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
RESULT

CHECK((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
	FeatureFinderAlgorithm<Peak1D,Feature>* ptr2 = FFASS::create();
	TEST_NOT_EQUAL(ptr2,0)
	delete ptr2;
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(FFASS::getProductName(),"simple")
RESULT

//remove log file
File::remove("featurefinder.log");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
