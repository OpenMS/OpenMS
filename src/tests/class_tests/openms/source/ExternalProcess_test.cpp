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
#include <OpenMS/SYSTEM/ExternalProcess.h>
///////////////////////////

#include <OpenMS/config.h>

using namespace OpenMS;
using namespace std;

// we just need ANY commandline tool available on (hopefully) all boxes.
// note that commands like "dir" or "type" are only known within cmd.exe and are not actual executables (unlike on Linux)
#ifdef OPENMS_WINDOWSPLATFORM
  const String exe = "cmd";
  const std::vector<String> args = {"/C", "echo hi"};
  const std::vector<String> args_broken = {"/C", "doesnotexist"};
#else
  const String exe = "ls";
  const std::vector<String> args = {"-l"};
  const std::vector<String> args_broken = {"-0"};
#endif //

START_TEST(ExternalProcess, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION(ExternalProcess())
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(~ExternalProcess())
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(void setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(RETURNSTATE run(const String& exe, const std::vector<String>& args, const String& working_dir, bool verbose, String& error_msg, IO_MODE io_mode, const std::map<String, String>& env, std::function<void()> idle_callback))
{
  String error_msg;
  { // without callbacks
    ExternalProcess ep;
    String error_msg;
    auto r = ep.run(exe, args, "", true, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::SUCCESS)
    TEST_EQUAL(error_msg.size(), 0)

    r = ep.run("this_exe_does_not_exist", args, "", true, error_msg);
    TEST_EQUAL(r,ExternalProcess::RETURNSTATE::FAILED_TO_START)
    TEST_NOT_EQUAL(error_msg.size(), 0);

    r = ep.run(exe, args_broken, "", true, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::NONZERO_EXIT)
    TEST_NOT_EQUAL(error_msg.size(), 0);
  }
  { // with callbacks
    String all_out, all_err;
    auto l_out = [&](const String& out) {all_out += out;};
    auto l_err = [&](const String& out) {all_err += out;};
    ExternalProcess ep(l_out, l_err);
    auto r = ep.run(exe, args, "", true, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::SUCCESS)
    TEST_EQUAL(error_msg.size(), 0);
    TEST_NOT_EQUAL(all_out.size(), 0)
    TEST_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();

    r = ep.run(exe, args_broken, "", false, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::NONZERO_EXIT)
    TEST_NOT_EQUAL(error_msg.size(), 0);
    TEST_EQUAL(all_out.size(), 0)
    std::cout << all_out << "\n\n";
    TEST_NOT_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();

    ep.setCallbacks(l_err, l_out); // swap callbacks
    r = ep.run(exe, args_broken, "", false, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::NONZERO_EXIT)
    TEST_NOT_EQUAL(error_msg.size(), 0);
    TEST_NOT_EQUAL(all_out.size(), 0)
    TEST_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();
  }
}
END_SECTION

START_SECTION(RETURNSTATE run(const String& exe, const std::vector<String>& args, const String& working_dir, bool verbose, IO_MODE io_mode, const std::map<String, String>& env, std::function<void()> idle_callback))
 NOT_TESTABLE // tested above..
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


