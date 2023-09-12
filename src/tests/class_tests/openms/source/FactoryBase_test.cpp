// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/FactoryBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FactoryBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FactoryBase* ptr = nullptr;
FactoryBase* nullPointer = nullptr;
START_SECTION(FactoryBase())
{
  ptr = new FactoryBase();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FactoryBase())
{
  delete ptr;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

