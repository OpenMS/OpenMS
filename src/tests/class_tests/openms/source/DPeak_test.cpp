// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/DPeak.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DPeak<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DPeak<1>::Type* ptr1 = nullptr;
DPeak<1>::Type* nullPointer1 = nullptr;
START_SECTION(DPeak())
{
	ptr1 = new DPeak<1>::Type();
  TEST_NOT_EQUAL(ptr1, nullPointer1);
}
END_SECTION

START_SECTION(~DPeak())
{
	delete ptr1;
}
END_SECTION

DPeak<2>::Type* ptr2 = nullptr;
DPeak<2>::Type* nullPointer2 = nullptr;
START_SECTION([EXTRA]DPeak())
{
	ptr2 = new DPeak<2>::Type();
  TEST_NOT_EQUAL(ptr2, nullPointer2);
}
END_SECTION

START_SECTION([EXTRA]~DPeak())
{
	delete ptr2;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
