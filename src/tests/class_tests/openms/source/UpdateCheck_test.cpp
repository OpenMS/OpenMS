// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/SYSTEM/UpdateCheck.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <fstream>
#include <sstream>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(UpdateCheck, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static void run(const std::string& tool_name, const std::string& version, int debug_level)))
{
  // UpdateCheck stamps a per-tool ".ver" file into File::getOpenMSConfigDir(). When that directory
  // cannot be created (a read-only / unwritable config location), run() must log a warning and return
  // *before* issuing any network request: it must not throw or crash, and the tool must continue.
  //
  // On unix we force an uncreatable config dir by pointing XDG_CONFIG_HOME at a path whose parent is a
  // regular file. getOpenMSConfigDir() then returns "<file>/OpenMS", and create_directories() fails with
  // ENOTDIR for every user (including root), so the test does not depend on filesystem permissions.
#ifdef __unix__
  const char* xdg_backup = getenv("XDG_CONFIG_HOME");

  // create a real temporary *file* and use it as the (bogus) config-home directory
  std::string blocker = File::getTemporaryFile();
  {
    std::ofstream f(blocker.c_str());
    f << "not a directory";
  }
  TEST_EQUAL(File::exists(blocker), true)
  setenv("XDG_CONFIG_HOME", blocker.c_str(), 1);

  // capture warnings emitted on the failure path
  std::ostringstream captured_warn;
  OPENMS_LOG_WARN.insert(captured_warn);

  bool threw = false;
  try
  {
    UpdateCheck::run("UpdateCheck_test_tool", "1.0.0", 0);
  }
  catch (...)
  {
    threw = true;
  }

  OPENMS_LOG_WARN.remove(captured_warn);

  // restore the environment so later tests are unaffected
  if (xdg_backup) { setenv("XDG_CONFIG_HOME", xdg_backup, 1); }
  else { unsetenv("XDG_CONFIG_HOME"); }

  // graceful degradation: no exception, and the documented warning was emitted
  TEST_EQUAL(threw, false)
  TEST_EQUAL(captured_warn.str().find("Could not create config directory") != std::string::npos, true)
#else
  NOT_TESTABLE // XDG_CONFIG_HOME redirection of the config dir is unix-specific
#endif
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
