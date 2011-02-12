// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmQT, "$Id FeatureGroupingAlgorithmQT_test.C 7941 2011-02-11 22:31:39Z hendrikweisser $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmQT* ptr = 0;
START_SECTION((FeatureGroupingAlgorithmQT()))
	ptr = new FeatureGroupingAlgorithmQT();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmQT()))
	delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
	FeatureGroupingAlgorithm* ptr2 = 0;
	ptr2 = FeatureGroupingAlgorithmQT::create();
	TEST_NOT_EQUAL(ptr2, 0)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(FeatureGroupingAlgorithmQT::getProductName(), "unlabeled_qt")
END_SECTION

START_SECTION((virtual void group(const vector< FeatureMap<> >&, ConsensusMap&)))
	// This is tested extensively in TEST/TOPP. See FeatureLinker_test.
	NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void group(const vector<ConsensusMap>&, ConsensusMap&)))
	// This is tested extensively in TEST/TOPP. See FeatureLinker_test.
	NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



