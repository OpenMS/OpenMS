// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/MzPeakFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MzPeakFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzPeakFile* ptr = nullptr;
MzPeakFile* nullPointer = nullptr;
START_SECTION((MzPeakFile()))
ptr = new MzPeakFile;
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzPeakFile()))
delete ptr;
END_SECTION

START_SECTION(void load(const String& filename, MapType& map))
{
  MSExperiment exp;
  // The handler is currently a stub and must throw NotImplemented.
  TEST_EXCEPTION(Exception::NotImplemented, MzPeakFile().load("nonexistent.mzpeak", exp))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
