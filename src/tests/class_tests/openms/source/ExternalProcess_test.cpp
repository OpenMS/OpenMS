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

#include <filesystem>
#include <fstream>

using namespace OpenMS;
using namespace std;

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

START_TEST(ExternalProcess, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION(ExternalProcess())
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(ExternalProcess(std::function<void(const std::string&)> callbackStdOut, std::function<void(const std::string&)> callbackStdErr))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(~ExternalProcess())
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(void setCallbacks(std::function<void(const std::string&)> callbackStdOut, std::function<void(const std::string&)> callbackStdErr))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(RETURNSTATE run(const std::string& exe, const std::vector<std::string>& args, const std::string& working_dir, bool verbose, std::string& error_msg, IO_MODE io_mode, const std::map<std::string, std::string>& env, std::function<void()> idle_callback))
{
  std::string error_msg;
  { // without callbacks
    ExternalProcess ep;
    std::string error_msg;
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
    std::string all_out, all_err;
    auto l_out = [&](const std::string& out) {all_out += out;};
    auto l_err = [&](const std::string& out) {all_err += out;};
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

START_SECTION(RETURNSTATE run(const std::string& exe, const std::vector<std::string>& args, const std::string& working_dir, bool verbose, IO_MODE io_mode, const std::map<std::string, std::string>& env, std::function<void()> idle_callback))
 NOT_TESTABLE // tested above..
END_SECTION

START_SECTION([EXTRA] run with spaces in the executable path and arguments)
{
#ifndef OPENMS_WINDOWSPLATFORM
  // An executable whose path contains spaces (e.g. "/opt/My Tool/bin/x") must
  // launch correctly and its stdout must be captured intact. Likewise a single
  // argument that itself contains spaces must be passed as ONE argument, not
  // shell-split. boost::process receives argv as a vector, so no shell re-quoting
  // should occur. (Windows is guarded out here because launching a freshly-written
  // .bat through ExternalProcess needs a cmd shim; the cmd-based section above
  // already exercises the Windows path.)
  std::string tmp;
  NEW_TMP_FILE(tmp)
  const std::filesystem::path dir = std::filesystem::path(tmp).parent_path() / "open ms space dir";
  std::filesystem::create_directories(dir);
  const std::filesystem::path script = dir / "my script.sh";
  {
    std::ofstream os(script);
    os << "#!/bin/sh\n"
          "echo MARKER_STDOUT_OK\n"
          "echo \"arg=[$1]\"\n";
  }
  std::filesystem::permissions(script, std::filesystem::perms::owner_all, std::filesystem::perm_options::add);

  std::string all_out, all_err, error_msg;
  ExternalProcess ep([&](const std::string& s) { all_out += s; },
                     [&](const std::string& s) { all_err += s; });

  // single argument that itself contains spaces
  const std::vector<std::string> spaced_args{"one two three"};
  auto r = ep.run(script.string(), spaced_args, "", true, error_msg);

  TEST_EQUAL(r, ExternalProcess::RETURNSTATE::SUCCESS)
  // stdout from the spaces-in-path executable was captured
  TEST_TRUE(all_out.find("MARKER_STDOUT_OK") != std::string::npos)
  // the spaced argument arrived as a single, intact argument (not split into 3)
  TEST_TRUE(all_out.find("arg=[one two three]") != std::string::npos)

  std::filesystem::remove_all(dir);
#else
  NOT_TESTABLE // Windows path-with-spaces is exercised via the cmd-based section above
#endif
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


