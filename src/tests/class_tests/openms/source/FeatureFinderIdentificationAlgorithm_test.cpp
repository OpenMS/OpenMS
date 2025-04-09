// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderIdentificationAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderIdentificationAlgorithm* ptr = 0;
FeatureFinderIdentificationAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderIdentificationAlgorithm())
{
  ptr = new FeatureFinderIdentificationAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFinderIdentificationAlgorithm())
{
  delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



