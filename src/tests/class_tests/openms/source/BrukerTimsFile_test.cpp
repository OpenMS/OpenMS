// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

using namespace OpenMS;
using namespace std;

START_TEST(BrukerTimsFile, "$Id$")

START_SECTION(void load(const String& path, MSExperiment& exp, const Config& config))
{
  // Test: invalid path throws FileNotReadable
  BrukerTimsFile f;
  MSExperiment exp;
  TEST_EXCEPTION(Exception::FileNotReadable, f.load("/nonexistent/path.d", exp));
}
END_SECTION

START_SECTION([FileHandler] BRUKER_TDF detection)
{
  // Test: .d suffix is detected as BRUKER_TDF
  TEST_EQUAL(FileHandler::getTypeByFileName("sample.d"), FileTypes::BRUKER_TDF);
  TEST_EQUAL(FileHandler::getTypeByFileName("/path/to/experiment.d"), FileTypes::BRUKER_TDF);

  // Test: non-.d suffixes are not BRUKER_TDF
  TEST_NOT_EQUAL(FileHandler::getTypeByFileName("sample.mzML"), FileTypes::BRUKER_TDF);
}
END_SECTION

END_TEST

#else // WITH_TIMSRUST

// Minimal test when timsrust is not available
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

START_TEST(BrukerTimsFile, "$Id$")
// No tests when WITH_TIMSRUST is off
END_TEST

#endif // WITH_TIMSRUST
