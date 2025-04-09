// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/INTERFACES/DataStructures.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DataStructures, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DataStructures* ptr = 0;
DataStructures* null_ptr = 0;
START_SECTION(DataStructures())
{
	ptr = new DataStructures();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~DataStructures())
{
	delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



