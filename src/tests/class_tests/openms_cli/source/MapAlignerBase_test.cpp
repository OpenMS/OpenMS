// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MapAlignerBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MapAlignerBase* ptr = 0;
MapAlignerBase* null_ptr = 0;
START_SECTION(MapAlignerBase())
{
	ptr = new MapAlignerBase();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MapAlignerBase())
{
	delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



