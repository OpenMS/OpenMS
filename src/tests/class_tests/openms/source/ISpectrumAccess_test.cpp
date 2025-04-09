// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/INTERFACES/ISpectrumAccess.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ISpectrumAccess, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ISpectrumAccess* ptr = 0;
ISpectrumAccess* null_ptr = 0;
START_SECTION(ISpectrumAccess())
{
	ptr = new ISpectrumAccess();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ISpectrumAccess())
{
	delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



