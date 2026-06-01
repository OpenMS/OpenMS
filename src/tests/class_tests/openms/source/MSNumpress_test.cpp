// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSNumpress, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSNumpress* ptr = 0;
MSNumpress* null_ptr = 0;
START_SECTION(MSNumpress())
{
	ptr = new MSNumpress();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MSNumpress())
{
	delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



