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

#include <cstdlib>
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
    : TOPPExternalToolBase("TOPPExternalToolBaseTest", "A test class", false, {}, false)
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
  const std::vector<std::string> args = {"-l"};
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

/////////////////////////////////////////////////////////////
END_TEST
