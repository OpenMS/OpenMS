// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderDefs, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderDefs* ptr = nullptr;
FeatureFinderDefs* null_ptr = nullptr;
START_SECTION(FeatureFinderDefs())
{
	ptr = new FeatureFinderDefs();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFinderDefs())
{
	delete ptr;
}
END_SECTION

START_SECTION(([FeatureFinderDefs::NoSuccessor] NoSuccessor(const char *file, int line, const char *function, const IndexPair &index)))
{
  // TODO
}
END_SECTION

START_SECTION(([FeatureFinderDefs::NoSuccessor] virtual ~NoSuccessor()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



