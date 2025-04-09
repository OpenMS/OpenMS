// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GlobalExceptionHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GlobalExceptionHandler* ptr = 0;
GlobalExceptionHandler* null_ptr = 0;
START_SECTION(GlobalExceptionHandler())
{
	ptr = new GlobalExceptionHandler();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~GlobalExceptionHandler())
{
	delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



