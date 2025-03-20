// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Joshua Charkow$
// $Authors: Joshua Charkow$
// --------------------------------------------------------------------------

#include "OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h"
#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SwathMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(testIsEqual)
{
	// map 1 and map 2 are equal
	OpenSwath::SwathMap map1;
	OpenSwath::SwathMap map2;
	TEST_EQUAL(map1.isEqual(map2), true);

	// map 3 and map 4 are equal 
	// map 5 is different because of ms1
	// map 6,7 is different because of mz bounds
	OpenSwath::SwathMap map3(1.0, 2.0, 1.5, false);
	OpenSwath::SwathMap map4(1.0, 2.0, 1.5, false);
	OpenSwath::SwathMap map5(1.0, 2.0, 1.5, true);
	TEST_EQUAL(map3.isEqual(map4), true);
	TEST_EQUAL(map3.isEqual(map5), false);

	// map 6,7 are different from map 3 different because of mz bounds
	OpenSwath::SwathMap map6(1.0, 3.0, 2.0, false);
	OpenSwath::SwathMap map7(2.0, 3.0, 2.5, false);

  // map 8 should be the same as map 3
	OpenSwath::SwathMap map8(1.0, 2.0, 1.5, -1, -1, false);
	TEST_EQUAL(map3.isEqual(map8), true);

  // map 9, 10 are equal
	OpenSwath::SwathMap map9(1.0, 2.0, 1.5, 1.0, 1.1, false);
	OpenSwath::SwathMap map10(1.0, 2.0, 1.5, 1.0, 1.1, false);
	TEST_EQUAL(map9.isEqual(map10), true);

  // map 11/12 is different from map 9 because of im bounds
	OpenSwath::SwathMap map11(1.0, 2.0, 1.5, 1.3, 1.4, false);
	OpenSwath::SwathMap map12(1.0, 2.0, 1.5, 1.0, 1.2, false);
	TEST_EQUAL(map9.isEqual(map11), false);
	TEST_EQUAL(map9.isEqual(map12), false);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST