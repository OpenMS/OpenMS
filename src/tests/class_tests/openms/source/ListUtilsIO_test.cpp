// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

using namespace OpenMS;
using namespace std;

START_TEST(ListUtilsIO, "$Id$")

START_SECTION(([EXTRA] template<typename StringType> StringList& operator<<(StringList& sl, const StringType& string)))
  StringList list;
  list << "a" << "b" << "c" << "a";
  TEST_EQUAL(list.size(),4)
  ABORT_IF(list.size() != 4)
  TEST_STRING_EQUAL(list[0],"a")
  TEST_STRING_EQUAL(list[1],"b")
  TEST_STRING_EQUAL(list[2],"c")
  TEST_STRING_EQUAL(list[3],"a")
END_SECTION

END_TEST
