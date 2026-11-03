// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/TOPPExternalToolBase.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

// minimal concrete tool that exposes the protected runExternalProcess_ for testing
class TOPPExternalToolBaseTest
  : public TOPPExternalToolBase
{
public:
  TOPPExternalToolBaseTest()
    : TOPPExternalToolBase("TOPPExternalToolBaseTest", "A test class", {}, false)
  {
    char* var = (char*)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
    _putenv(var);
#else
    putenv(var);
#endif
    main(0, nullptr);
  }

  void registerOptionsAndFlags_() override
  {
  }

  ExitCodes main_(int /*argc*/, const char** /*argv*/) override
  {
    return EXECUTION_OK;
  }

  TOPPBase::ExitCodes runExternalProcess(const std::string& executable, const std::vector<std::string>& arguments, const std::string& workdir) const
  {
    return runExternalProcess_(executable, arguments, workdir);
  }

  bool findExecutable(std::string& executable) const
  {
    return findExecutable_(executable);
  }
};

START_TEST(TOPPExternalToolBase, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] ExitCodes runExternalProcess_(const std::string& executable, const std::vector<std::string>& arguments, const std::string& workdir) const))
{

// we just need ANY commandline tool available on (hopefully) all boxes.
// note that commands like "dir" or "type" are only known within cmd.exe and are not actual executables (unlike on Linux)
#ifdef OPENMS_WINDOWSPLATFORM
  const std::string exe = "cmd";
  const std::vector<std::string> args = {"/C", "echo hi"};
  const std::vector<std::string> args_broken = {"/C", "doesnotexist"};
#else
  const std::string exe = "ls";
  // Inspect the directory itself: parallel tests may delete entries while ls lists them.
  const std::vector<std::string> args = {"-ld", "."};
  const std::vector<std::string> args_broken = {"-0"};
#endif //

  TOPPExternalToolBaseTest topp;
  auto result = topp.runExternalProcess("/path/does/not/exists.exe", {}, "");
  TEST_EQUAL(result, TOPPBase::EXTERNAL_PROGRAM_NOTFOUND);

  result = topp.runExternalProcess(exe, args_broken, "");
  TEST_EQUAL(result, TOPPBase::EXTERNAL_PROGRAM_ERROR);

  result = topp.runExternalProcess(exe, args, "");
  TEST_EQUAL(result, TOPPBase::EXECUTION_OK);
}
END_SECTION

START_SECTION((static StringList getThirdPartyToolLocations(const std::string& data_path)))
{
  TempDir tdir;
  std::string base = tdir.getPath();
  File::makeDir(base + "/THIRDPARTY/Sage");
  File::makeDir(base + "/THIRDPARTY/Comet");
  { // a file, not a tool folder
    std::ofstream f(std::string(base + "/THIRDPARTY/README.txt"));
    f << "test";
  }
  StringList tools = TOPPExternalToolBase::getThirdPartyToolLocations(base);
  TEST_EQUAL(tools.size(), 2)
  // sorted, with '/' as separator and at the end
  TEST_TRUE(StringUtils::hasSuffix(tools[0], "/THIRDPARTY/Comet/"))
  TEST_TRUE(StringUtils::hasSuffix(tools[1], "/THIRDPARTY/Sage/"))

  // no THIRDPARTY folder
  TEST_EQUAL(TOPPExternalToolBase::getThirdPartyToolLocations(base + "/THIRDPARTY/Sage").size(), 0)

  // a THIRDPARTY folder with only a file in it, as in a build tree
  File::makeDir(base + "/share/THIRDPARTY");
  {
    std::ofstream f(std::string(base + "/share/THIRDPARTY/ReadMe.txt"));
    f << "test";
  }
  TEST_EQUAL(TOPPExternalToolBase::getThirdPartyToolLocations(base + "/share").size(), 0)
}
END_SECTION

START_SECTION((static StringList getThirdPartyToolLocations()))
{
  // those of the shared-data directory this build uses
  TEST_TRUE(TOPPExternalToolBase::getThirdPartyToolLocations() == TOPPExternalToolBase::getThirdPartyToolLocations(File::getOpenMSDataPath()))
}
END_SECTION

START_SECTION((static bool findThirdPartyExecutable(std::string& exe_filename, const StringList& locations)))
{
  TempDir tdir;
  std::string base = tdir.getPath();
  File::makeDir(base + "/THIRDPARTY/Comet");
  File::makeDir(base + "/THIRDPARTY/Sage/sage"); // a folder named like the executable
  {
    std::ofstream f(std::string(base + "/THIRDPARTY/Comet/comet.exe"));
    f << "test";
  }
  StringList locations = TOPPExternalToolBase::getThirdPartyToolLocations(base);
  TEST_EQUAL(locations.size(), 2)

  std::string exe = "comet.exe";
  TEST_TRUE(TOPPExternalToolBase::findThirdPartyExecutable(exe, locations))
  TEST_EQUAL(exe, locations[0] + "comet.exe")

  // a folder is not an executable
  exe = "sage";
  TEST_FALSE(TOPPExternalToolBase::findThirdPartyExecutable(exe, locations))
  TEST_EQUAL(exe, "sage")

  // not there: the name stays as it is
  exe = "percolator";
  TEST_FALSE(TOPPExternalToolBase::findThirdPartyExecutable(exe, locations))
  TEST_EQUAL(exe, "percolator")

  // a name with a directory part is not looked up, even where it would exist
  exe = "../Comet/comet.exe";
  TEST_FALSE(TOPPExternalToolBase::findThirdPartyExecutable(exe, locations))
  TEST_EQUAL(exe, "../Comet/comet.exe")
  exe = "Comet\\comet.exe";
  TEST_FALSE(TOPPExternalToolBase::findThirdPartyExecutable(exe, StringList{base + "/THIRDPARTY/"}))

  exe = "";
  TEST_FALSE(TOPPExternalToolBase::findThirdPartyExecutable(exe, locations))

  // no folders
  exe = "comet.exe";
  TEST_FALSE(TOPPExternalToolBase::findThirdPartyExecutable(exe, StringList()))
  TEST_EQUAL(exe, "comet.exe")
}
END_SECTION

START_SECTION((static bool findThirdPartyExecutable(std::string& exe_filename)))
{
  // searches the folders of this build's shared-data directory
  std::string exe = "comet.exe", expected = "comet.exe";
  TEST_EQUAL(TOPPExternalToolBase::findThirdPartyExecutable(exe), TOPPExternalToolBase::findThirdPartyExecutable(expected, TOPPExternalToolBase::getThirdPartyToolLocations()))
  TEST_EQUAL(exe, expected)

  // a name with a directory part is left as it is
  exe = "/does/not/exist/comet.exe";
  TEST_FALSE(TOPPExternalToolBase::findThirdPartyExecutable(exe))
  TEST_EQUAL(exe, "/does/not/exist/comet.exe")
}
END_SECTION

START_SECTION(([EXTRA] bool findExecutable_(std::string& executable) const))
{
  TOPPExternalToolBaseTest topp;
#ifdef OPENMS_WINDOWSPLATFORM
  std::string exe = "cmd";
#else
  std::string exe = "ls";
#endif
  // on the PATH: resolved to its full path
  TEST_TRUE(topp.findExecutable(exe))
  TEST_TRUE(exe.find_first_of("/\\") != std::string::npos)
  TEST_TRUE(File::exists(exe))

  // neither on the PATH nor among the third-party tools
  exe = "does_not_exist_anywhere_4711";
  TEST_FALSE(topp.findExecutable(exe))
  TEST_EQUAL(exe, "does_not_exist_anywhere_4711")
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
