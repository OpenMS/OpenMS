// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SYSTEM/SystemSettings.h>
///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>

using namespace OpenMS;
using namespace std;

namespace
{
  // getenv() pointers may be invalidated by later setenv() calls, so copy the value
  std::optional<std::string> readEnv(const char* name)
  {
    const char* value = getenv(name);
    if (value == nullptr) return std::nullopt;
    return std::string(value);
  }

  void setEnv(const char* name, const std::string& value)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    _putenv_s(name, value.c_str());
#else
    setenv(name, value.c_str(), 1);
#endif
  }

  void restoreEnv(const char* name, const std::optional<std::string>& previous)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    _putenv_s(name, previous ? previous->c_str() : "");
#else
    if (previous) setenv(name, previous->c_str(), 1);
    else unsetenv(name);
#endif
  }
}

START_TEST(SystemSettings, "$Id$")

// remembered before any section changes it; restored in the getOpenMSHomePath() section
const std::optional<std::string> home_backup = readEnv("OPENMS_HOME_PATH");
// the directory OPENMS_HOME_PATH is pointed at by the getUserDirectory() section
std::string home_dirname;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::string getTempDirectory())
  TEST_NOT_EQUAL(SystemSettings::getTempDirectory(), std::string())
  TEST_EQUAL(File::exists(SystemSettings::getTempDirectory()), true)
END_SECTION

START_SECTION(static std::string getUserDirectory())
  TEST_NOT_EQUAL(SystemSettings::getUserDirectory(), std::string())
  TEST_EQUAL(File::exists(SystemSettings::getUserDirectory()), true)

  // Point OPENMS_HOME_PATH at a fresh directory and check that it takes precedence
  // (no changes on the file system occur). The previous value is restored at the
  // end of the getOpenMSHomePath() section below.
  home_dirname = SystemSettings::getTempDirectory() + "/" + File::getUniqueName() + "/";
  TEST_EQUAL(File::makeDir(home_dirname), true);
  setEnv("OPENMS_HOME_PATH", home_dirname);
  TEST_EQUAL(SystemSettings::getUserDirectory(), home_dirname)
  // Note: this does not guarantee any more that the user directory or an
  // OpenMS.ini file exists at the new location.
END_SECTION

START_SECTION(static std::string getOpenMSHomePath())
  // OPENMS_HOME_PATH was set in the section above and takes precedence
  TEST_NOT_EQUAL(SystemSettings::getOpenMSHomePath(), std::string())
  TEST_EQUAL(SystemSettings::getOpenMSHomePath(), home_dirname)
  restoreEnv("OPENMS_HOME_PATH", home_backup);
END_SECTION

START_SECTION(static std::string getOpenMSConfigDir())
  std::string config_dir = SystemSettings::getOpenMSConfigDir();
  TEST_NOT_EQUAL(config_dir, std::string())
  // every platform branch resolves to a folder named "OpenMS" with no trailing separator
  TEST_EQUAL(StringUtils::hasSuffix(config_dir, "OpenMS"), true)
  TEST_EQUAL(StringUtils::hasSuffix(config_dir, "/"), false)
#ifdef __unix__
  // on unix-like systems, XDG_CONFIG_HOME takes precedence when set
  const std::optional<std::string> xdg_backup = readEnv("XDG_CONFIG_HOME");
  setEnv("XDG_CONFIG_HOME", "/tmp/openms_xdg_test");
  TEST_EQUAL(SystemSettings::getOpenMSConfigDir(), "/tmp/openms_xdg_test/OpenMS")
  // restore previous environment to avoid side effects on later tests
  restoreEnv("XDG_CONFIG_HOME", xdg_backup);
#endif
END_SECTION

START_SECTION(static Param getSystemParameters())
  Param p = SystemSettings::getSystemParameters();
  TEST_EQUAL(!p.empty(), true)
  TEST_EQUAL(p.getValue("version"), VersionInfo::getVersion())
  // the defaults always carry the directory policy entries, even without an OpenMS.ini
  TEST_TRUE(p.exists("home_dir"))
  TEST_TRUE(p.exists("temp_dir"))
  TEST_TRUE(p.exists("id_db_dir"))
END_SECTION

START_SECTION(static std::string findDatabase(const std::string& db_name))
  // findDatabase() logs the miss before rethrowing; the exception is what this asserts on.
  {
    Logger::LogSinkGuard quiet(getThreadLocalLogError(), std::cerr);
    TEST_EXCEPTION(Exception::FileNotFound, SystemSettings::findDatabase("filedoesnotexists"))
  }
  // The success path is chatty too -- it announces the resolved path on the info log (which goes
  // to stdout, so the stderr check above does not cover it), and since nothing flushes it until
  // teardown it surfaces *after* the test's PASSED banner. Same packaging-log noise as the errors.
  std::string db;
  {
    Logger::LogSinkGuard quiet(getThreadLocalLogInfo(), std::cout);
    db = SystemSettings::findDatabase("./CV/unimod.obo");
  }
  TEST_EQUAL(StringUtils::hasSubstring(db, "share/OpenMS"), true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
