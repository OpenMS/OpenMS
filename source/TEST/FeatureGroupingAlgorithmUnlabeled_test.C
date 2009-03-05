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
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmUnlabeled, "$Id FeatureFinder_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmUnlabeled* ptr = 0;
START_SECTION((FeatureGroupingAlgorithmUnlabeled()))
	ptr = new FeatureGroupingAlgorithmUnlabeled();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmUnlabeled()))
	delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
	FeatureGroupingAlgorithm* ptr2 = 0;
	ptr2 = FeatureGroupingAlgorithmUnlabeled::create();
	TEST_NOT_EQUAL(ptr2, 0)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(FeatureGroupingAlgorithmUnlabeled::getProductName(),"unlabeled")
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap<> > &maps, ConsensusMap &out)))
	// This is tested extensively in TEST/TOPP. See FeatureLinker_test.
	NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



