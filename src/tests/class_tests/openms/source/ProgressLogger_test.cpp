// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/ProgressLogger.h>

///////////////////////////

START_TEST(ProgressLogger, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ProgressLogger* ptr = nullptr;
ProgressLogger* null_ptr = nullptr;

START_SECTION((ProgressLogger()))
{
  ptr = new ProgressLogger();
  TEST_NOT_EQUAL(ptr, null_ptr)
  // the default log type is NONE (progress logging disabled)
  TEST_EQUAL(ptr->getLogType(), ProgressLogger::NONE)
}
END_SECTION

START_SECTION((virtual ~ProgressLogger()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void setLogType(LogType type) const))
{
  ProgressLogger pl;
  pl.setLogType(ProgressLogger::CMD);
  TEST_EQUAL(pl.getLogType(), ProgressLogger::CMD)
  pl.setLogType(ProgressLogger::GUI);
  TEST_EQUAL(pl.getLogType(), ProgressLogger::GUI)
  pl.setLogType(ProgressLogger::NONE);
  TEST_EQUAL(pl.getLogType(), ProgressLogger::NONE)
}
END_SECTION

START_SECTION((LogType getLogType() const))
{
  ProgressLogger pl;
  TEST_EQUAL(pl.getLogType(), ProgressLogger::NONE)
  pl.setLogType(ProgressLogger::CMD);
  TEST_EQUAL(pl.getLogType(), ProgressLogger::CMD)
}
END_SECTION

START_SECTION((ProgressLogger(const ProgressLogger& other)))
{
  ProgressLogger pl;
  pl.setLogType(ProgressLogger::CMD);
  ProgressLogger copy(pl);
  TEST_EQUAL(copy.getLogType(), ProgressLogger::CMD)

  ProgressLogger pl2;
  pl2.setLogType(ProgressLogger::GUI);
  ProgressLogger copy2(pl2);
  TEST_EQUAL(copy2.getLogType(), ProgressLogger::GUI)
}
END_SECTION

START_SECTION((ProgressLogger & operator=(const ProgressLogger& other)))
{
  ProgressLogger pl;
  pl.setLogType(ProgressLogger::CMD);
  ProgressLogger assigned;
  assigned = pl;
  TEST_EQUAL(assigned.getLogType(), ProgressLogger::CMD)
}
END_SECTION

START_SECTION((void startProgress(SignedSize begin, SignedSize end, const std::string& label) const))
{
  // With logging disabled (NONE) the progress calls are safe no-ops.
  ProgressLogger pl;
  pl.setLogType(ProgressLogger::NONE);
  pl.startProgress(0, 10, "test");
  pl.setProgress(5);
  pl.endProgress();
  TEST_EQUAL(pl.getLogType(), ProgressLogger::NONE)
}
END_SECTION

START_SECTION((void setProgress(SignedSize value) const))
{
  ProgressLogger pl;
  pl.setLogType(ProgressLogger::NONE);
  pl.startProgress(0, 100, "range");
  pl.setProgress(0);
  pl.setProgress(50);
  pl.setProgress(100);
  pl.endProgress();
  NOT_TESTABLE // no observable state in NONE mode; verifies the calls do not crash
}
END_SECTION

START_SECTION((void nextProgress() const))
{
  ProgressLogger pl;
  pl.setLogType(ProgressLogger::NONE);
  pl.startProgress(0, 3, "next");
  pl.nextProgress();
  pl.nextProgress();
  pl.endProgress();
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void endProgress(UInt64 bytes_processed = 0) const))
{
  ProgressLogger pl;
  pl.setLogType(ProgressLogger::NONE);
  pl.startProgress(0, 10, "bytes");
  pl.endProgress(1024); // optional byte count for MB/s estimate (ignored in NONE mode)
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
