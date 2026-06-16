// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/RNaseDB.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>

#include <algorithm>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(RNaseDB, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static const RNaseDB* getInstance()))
{
  const RNaseDB* db = RNaseDB::getInstance();
  TEST_NOT_EQUAL(db, nullptr)
  // singleton: repeated calls return the same instance
  TEST_EQUAL(db == RNaseDB::getInstance(), true)
}
END_SECTION

START_SECTION([EXTRA] hasEnzyme / getEnzyme (inherited from DigestionEnzymeDB))
{
  const RNaseDB* db = RNaseDB::getInstance();
  // known RNases
  TEST_EQUAL(db->hasEnzyme("RNase_T1"), true)
  TEST_EQUAL(db->hasEnzyme("RNase_U2"), true)
  TEST_EQUAL(db->hasEnzyme("cusativin"), true)
  // not a known RNase
  TEST_EQUAL(db->hasEnzyme("NotAnRNase"), false)

  const DigestionEnzymeRNA* e = db->getEnzyme("RNase_T1");
  TEST_NOT_EQUAL(e, nullptr)
  TEST_STRING_EQUAL(e->getName(), "RNase_T1")

  // an unknown enzyme throws
  TEST_EXCEPTION(Exception::ElementNotFound, db->getEnzyme("NotAnRNase"))
}
END_SECTION

START_SECTION([EXTRA] getAllNames (inherited from DigestionEnzymeDB))
{
  const RNaseDB* db = RNaseDB::getInstance();
  std::vector<std::string> names;
  db->getAllNames(names);
  TEST_EQUAL(names.empty(), false)
  TEST_EQUAL(std::find(names.begin(), names.end(), std::string("RNase_T1")) != names.end(), true)
  TEST_EQUAL(std::find(names.begin(), names.end(), std::string("cusativin")) != names.end(), true)
  // every reported name is recognized by hasEnzyme
  for (const auto& n : names)
  {
    TEST_EQUAL(db->hasEnzyme(n), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
