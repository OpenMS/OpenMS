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
#include <OpenMS/SYSTEM/RWrapper.h>
///////////////////////////

#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(RWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// RWrapper drives the 'Rscript' interpreter through boost::process. These tests
// pin the graceful-degradation contract that callers rely on: missing interpreter
// or missing script must be reported as 'false' / FileNotFound, never an uncaught
// exception or crash. The paths below are deterministic and need no R installed
// (the "skip if Rscript absent" behavior). The R-present positive path is guarded
// by findR() and only runs where 'Rscript' is on PATH.

START_SECTION((static bool findR(const std::string& executable = "Rscript", bool verbose = true)))
{
  // A non-existent interpreter is reported as "not found" cleanly: boost::process
  // throws process_error internally, which findR catches and turns into 'false'.
  bool threw = false;
  bool found_bogus = true;
  try { found_bogus = RWrapper::findR("this_is_not_a_real_R_interpreter_xyz", false); }
  catch (...) { threw = true; }
  TEST_EQUAL(threw, false)
  TEST_EQUAL(found_bogus, false)

  // Probing for the real 'Rscript' returns a clean bool either way -- this is the
  // "is R available?" signal callers use to skip optional R-based steps.
  threw = false;
  try { (void)RWrapper::findR("Rscript", false); }
  catch (...) { threw = true; }
  TEST_EQUAL(threw, false)
}
END_SECTION

START_SECTION((static std::string findScript(const std::string& script_file, bool verbose = true)))
{
  // A non-existent script raises a clear FileNotFound -- this is the only RWrapper
  // method that throws (runScript() swallows it into a 'false' return; see below).
  TEST_EXCEPTION(Exception::FileNotFound,
                 RWrapper::findScript("definitely_nonexistent_script_qwerty.R", false))
}
END_SECTION

START_SECTION((static bool runScript(const std::string& script_file, const std::vector<std::string>& cmd_args, const std::string& executable = "Rscript", bool find_R = false, bool verbose = true)))
{
  const std::vector<std::string> no_args;

  // R absent (find_R=true with a bogus interpreter): runScript returns false cleanly.
  bool threw = false;
  bool ok = true;
  try { ok = RWrapper::runScript("any.R", no_args, "this_is_not_a_real_R_interpreter_xyz", true, false); }
  catch (...) { threw = true; }
  TEST_EQUAL(threw, false)
  TEST_EQUAL(ok, false)

  // Missing script (find_R=false): findScript throws FileNotFound internally, which
  // runScript catches and turns into a clean 'false' (no exception escapes).
  threw = false; ok = true;
  try { ok = RWrapper::runScript("definitely_nonexistent_script_qwerty.R", no_args, "Rscript", false, false); }
  catch (...) { threw = true; }
  TEST_EQUAL(threw, false)
  TEST_EQUAL(ok, false)

  // When 'Rscript' IS available, the find_R path resolves it, but a missing script
  // still fails gracefully (returns false). On systems without 'Rscript' this block
  // is skipped -- exactly the "skip if Rscript absent" behavior.
  if (RWrapper::findR("Rscript", false))
  {
    threw = false; ok = true;
    try { ok = RWrapper::runScript("definitely_nonexistent_script_qwerty.R", no_args, "Rscript", true, false); }
    catch (...) { threw = true; }
    TEST_EQUAL(threw, false)
    TEST_EQUAL(ok, false)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
