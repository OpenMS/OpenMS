// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/JavaInfo.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static bool canRun(const String &file)))
  // test for missing java executable
  TEST_EQUAL(JavaInfo::canRun(""), false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
