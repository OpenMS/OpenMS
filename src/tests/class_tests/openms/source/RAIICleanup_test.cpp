// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/RAIICleanup.h>

#include <functional>

///////////////////////////

START_TEST(RAIICleanup, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((RAIICleanup(std::function<void()> l)))
{
  // the cleanup lambda must not run before destruction
  bool ran = false;
  {
    RAIICleanup c([&]() { ran = true; });
    TEST_EQUAL(ran, false)
  } // c goes out of scope here -> lambda runs
  TEST_EQUAL(ran, true)
}
END_SECTION

START_SECTION((~RAIICleanup()))
{
  // lambda runs exactly once per instance, at end of scope
  int counter = 0;
  {
    RAIICleanup c([&]() { counter += 5; });
  }
  TEST_EQUAL(counter, 5)

  {
    RAIICleanup c([&]() { counter += 5; });
  }
  TEST_EQUAL(counter, 10)

  // the lambda can capture and mutate arbitrary state
  std::string s = "start";
  {
    RAIICleanup c([&]() { s = "cleaned"; });
    TEST_STRING_EQUAL(s, "start")
  }
  TEST_STRING_EQUAL(s, "cleaned")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
