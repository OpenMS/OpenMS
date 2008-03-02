// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marcel Grunert $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

//#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmWavelet.h>

///////////////////////////

START_TEST(FeatureFinderAlgorithmWavelet, "$Id: FeatureFinderAlgorithmWavelet_test.C 2783 2008-02-28 22:18:12Z grunert $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/*
typedef FeatureFinderAlgorithmWavelet<RawDataPoint1D,Feature> FFAW;

FFAW* ptr;
CHECK(FeatureFinderAlgorithmWavelet())
	ptr = new FFAW;
	TEST_NOT_EQUAL(ptr,0)
	delete ptr;
RESULT

CHECK(virtual void run())
RESULT

CHECK((static FeatureFinderAlgorithm<PeakType,FeatureType>* create()))
	FeatureFinderAlgorithm<RawDataPoint1D,Feature>* ptr2 = FFAW::create();
	TEST_NOT_EQUAL(ptr2,0)
	delete ptr2;
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(FFAW::getProductName(),"wavelet")
RESULT
*/
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
