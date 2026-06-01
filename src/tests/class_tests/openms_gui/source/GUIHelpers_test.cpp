// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GUIHelpers, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION(([EXTRA] size_t OverlapDetector::placeItem(double x_start, double x_end)))
	GUIHelpers::OverlapDetector od(3);
	TEST_EQUAL(od.placeItem(1, 3), 0);
  TEST_EQUAL(od.placeItem(1, 3), 1);
  TEST_EQUAL(od.placeItem(1, 3), 2);
  TEST_EQUAL(od.placeItem(1, 3), 0);

  TEST_EQUAL(od.placeItem(4, 8), 0);
  TEST_EQUAL(od.placeItem(5, 11), 1);
  TEST_EQUAL(od.placeItem(9, 10), 0);

  TEST_EQUAL(od.placeItem(12, 20), 0);
  TEST_EQUAL(od.placeItem(12, 18), 1);
  TEST_EQUAL(od.placeItem(12, 19), 2);

  TEST_EQUAL(od.placeItem(16, 25), 1);
  TEST_EQUAL(od.placeItem(16, 25), 2);
  TEST_EQUAL(od.placeItem(16, 25), 0);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



