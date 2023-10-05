// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/Utils/MatrixUtils.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MatrixUtils, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MatrixUtils* ptr = 0;
MatrixUtils* null_ptr = 0;
START_SECTION(MatrixUtils())
{
	ptr = new MatrixUtils();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MatrixUtils())
{
	delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



