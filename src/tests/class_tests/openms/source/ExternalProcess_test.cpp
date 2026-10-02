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

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <OpenMS/config.h>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <thread>

using namespace OpenMS;
using namespace std;

namespace
{
std::string self_executable;

std::string pathString(const std::filesystem::path& path)
{
  const auto utf8 = path.u8string();
  return std::string(utf8.begin(), utf8.end());
}
} // namespace

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

// Keep the native executable as a child-process fixture. In particular, a shell
// has different Windows quoting rules and cannot verify an argv round trip.
#define main externalProcessTestsMain
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
  // run everything in a private working directory: the test binary runs in a directory
  // shared by all tests, and under 'ctest --parallel' the other tests' temp files appear
  // and vanish while 'ls -l' walks it, making it print errors and exit non-zero (#9948)
  TempDir tmp_dir;
  { // one stable entry, so 'ls -l' has something to list
    std::ofstream file(tmp_dir.getPath() + "some_file.txt");
    file << "content\n";
  }

  std::string error_msg;
  { // without callbacks
    ExternalProcess ep;
    std::string error_msg;
    auto r = ep.run(exe, args, tmp_dir.getPath(), true, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::SUCCESS)
    TEST_EQUAL(error_msg.size(), 0)

    r = ep.run("this_exe_does_not_exist", args, tmp_dir.getPath(), true, error_msg);
    TEST_EQUAL(r,ExternalProcess::RETURNSTATE::FAILED_TO_START)
    TEST_NOT_EQUAL(error_msg.size(), 0);

    r = ep.run(exe, args_broken, tmp_dir.getPath(), true, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::NONZERO_EXIT)
    TEST_NOT_EQUAL(error_msg.size(), 0);
  }
  { // with callbacks
    std::string all_out, all_err;
    auto l_out = [&](const std::string& out) {all_out += out;};
    auto l_err = [&](const std::string& out) {all_err += out;};
    ExternalProcess ep(l_out, l_err);
    auto r = ep.run(exe, args, tmp_dir.getPath(), true, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::SUCCESS)
    TEST_EQUAL(error_msg.size(), 0);
    TEST_NOT_EQUAL(all_out.size(), 0)
    TEST_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();

    r = ep.run(exe, args_broken, tmp_dir.getPath(), false, error_msg);
    TEST_EQUAL(r, ExternalProcess::RETURNSTATE::NONZERO_EXIT)
    TEST_NOT_EQUAL(error_msg.size(), 0);
    TEST_EQUAL(all_out.size(), 0)
    std::cout << all_out << "\n\n";
    TEST_NOT_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();

    ep.setCallbacks(l_err, l_out); // swap callbacks
    r = ep.run(exe, args_broken, tmp_dir.getPath(), false, error_msg);
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

START_SECTION([EXTRA] native argument round trip preserves empty arguments quotes tabs and trailing backslashes)
{
  const std::vector<std::string> values {"",       "plain",      "one two",     "one\ttwo", "C:\\Data Set\\",       "C:\\Data Set\\\\", "a\"b",
                                         "a\\\"b", "trailing\\", "after empty", "",         "\xE6\xB8\xAC \xC3\xA4"};
  std::vector<std::string> arguments {"--openms-process-argv"};
  arguments.insert(arguments.end(), values.begin(), values.end());
  std::string output;
  ExternalProcess process([&](const std::string& text) { output += text; }, [](const std::string&) {});
  const auto result = process.runWithResult(self_executable, arguments, "", false);
  std::string expected;
  for (const auto& value : values)
    expected += std::to_string(value.size()) + ":" + value;
  TEST_EQUAL(result.state, ExternalProcess::RETURNSTATE::SUCCESS)
  TEST_EQUAL(output, expected)
}
END_SECTION

START_SECTION([EXTRA] structured results preserve numeric exit status and launch errors)
{
  ExternalProcess ep;
#ifdef OPENMS_WINDOWSPLATFORM
  const std::string shell = "cmd";
  const std::vector<std::string> exit_args = {"/C", "exit /B 37"};
#else
  const std::string shell = "/bin/sh";
  const std::vector<std::string> exit_args = {"-c", "exit 37"};
#endif
  for (const auto mode : {ExternalProcess::IO_MODE::READ_WRITE, ExternalProcess::IO_MODE::NO_IO})
  {
    const auto result = ep.runWithResult(shell, exit_args, "", false, mode);
    TEST_EQUAL(result.state, ExternalProcess::RETURNSTATE::NONZERO_EXIT)
    TEST_EQUAL(result.exit_code, 37)
    TEST_TRUE(result.error_message.find("37") != std::string::npos)
  }
  const auto missing = ep.runWithResult("openms_external_process_executable_that_does_not_exist", {}, "", false);
  TEST_EQUAL(missing.state, ExternalProcess::RETURNSTATE::FAILED_TO_START)
  TEST_EQUAL(missing.exit_code, -1)
  TEST_FALSE(missing.error_message.empty())

  std::atomic_bool cancel {true};
  const auto cancelled = ep.runWithResult(shell, exit_args, "", false, ExternalProcess::IO_MODE::READ_WRITE, {}, nullptr, &cancel);
  TEST_EQUAL(cancelled.state, ExternalProcess::RETURNSTATE::CANCELLED)
  TEST_EQUAL(cancelled.exit_code, -1)
}
END_SECTION

START_SECTION([EXTRA] Unicode executable arguments working directory and environment)
{
  TempDir tmp_dir;
  const auto directory = to_path(tmp_dir.getPath()) / to_path("unicode \xE6\xB8\xAC \xC3\xA4");
  std::filesystem::create_directories(directory);
  const std::string input_name = "input \xE6\xB8\xAC \xC3\xA4.txt";
  const std::string output_name = "copied \xE6\xB8\xAC \xC3\xA4.txt";
  {
    std::ofstream stream(directory / to_path(input_name));
    stream << "UNICODE_ARGUMENTS_OK";
  }
#ifdef OPENMS_WINDOWSPLATFORM
  const auto executable = directory / to_path("tool \xE6\xB8\xAC \xC3\xA4.exe");
  std::filesystem::copy_file(to_path(self_executable), executable);
  const std::vector<std::string> arguments {"--openms-process-copy-env", input_name};
#else
  const auto executable = directory / to_path("tool \xE6\xB8\xAC \xC3\xA4.sh");
  {
    std::ofstream stream(executable);
    stream << "#!/bin/sh\ncp -- \"$1\" \"$OPENMS_UTF8_OUTPUT\"\n";
  }
  std::filesystem::permissions(executable, std::filesystem::perms::owner_all, std::filesystem::perm_options::add);
  const std::vector<std::string> arguments {input_name};
#endif
  std::map<std::string, std::string> environment {{"OPENMS_UTF8_OUTPUT", output_name}};
#ifdef OPENMS_WINDOWSPLATFORM
  // The copied test executable must still find the DLLs beside the original.
  const char* inherited_path = std::getenv("PATH");
  environment["PATH"] = pathString(to_path(self_executable).parent_path()) + ";" + (inherited_path ? inherited_path : "");
#endif
  ExternalProcess process;
  const auto result
    = process.runWithResult(pathString(executable), arguments, pathString(directory), false, ExternalProcess::IO_MODE::READ_WRITE, environment);
  TEST_EQUAL(result.state, ExternalProcess::RETURNSTATE::SUCCESS)
  TEST_EQUAL(result.exit_code, 0)
  TEST_TRUE(std::filesystem::is_regular_file(directory / to_path(output_name)))
  std::ifstream stream(directory / to_path(output_name));
  std::string content;
  std::getline(stream, content);
  TEST_EQUAL(content, "UNICODE_ARGUMENTS_OK")
}
END_SECTION

START_SECTION([EXTRA] structured execution streams both channels and preserves environment and working directory)
{
  TempDir tmp_dir;
  std::string all_out, all_err;
  const char* inherited_value = std::getenv("OPENMS_EXTERNAL_PROCESS_TEST");
  const bool inherited_present = inherited_value != nullptr;
  const std::string original_value = inherited_present ? inherited_value : "";
  ExternalProcess ep([&](const std::string& text) { all_out += text; }, [&](const std::string& text) { all_err += text; });
#ifdef OPENMS_WINDOWSPLATFORM
  const std::string shell = "cmd";
  const std::vector<std::string> stream_args
    = {"/C",
       "echo %OPENMS_EXTERNAL_PROCESS_TEST% & echo STDERR_MARKER 1>&2 & echo CREATED>relative.txt & for /L %i in (1,1,10000) do @echo OUTPUT_LINE"};
#else
  const std::string shell = "/bin/sh";
  const std::vector<std::string> stream_args
    = {"-c", "printf '%s\\n' \"$OPENMS_EXTERNAL_PROCESS_TEST\"; printf 'STDERR_MARKER\\n' >&2; printf CREATED > relative.txt; i=0; while [ $i -lt "
             "10000 ]; do printf 'OUTPUT_LINE\\n'; i=$((i+1)); done"};
#endif
  const auto result = ep.runWithResult(shell, stream_args, tmp_dir.getPath(), false, ExternalProcess::IO_MODE::READ_WRITE,
                                       {{"OPENMS_EXTERNAL_PROCESS_TEST", "ENVIRONMENT_MARKER"}});
  TEST_EQUAL(result.state, ExternalProcess::RETURNSTATE::SUCCESS)
  TEST_EQUAL(result.exit_code, 0)
  TEST_TRUE(result.error_message.empty())
  TEST_TRUE(all_out.find("ENVIRONMENT_MARKER") != std::string::npos)
  TEST_TRUE(all_err.find("STDERR_MARKER") != std::string::npos)
  TEST_TRUE(std::filesystem::exists(std::filesystem::path(tmp_dir.getPath()) / "relative.txt"))
  // More than a pipe buffer must be drained while the child is still running.
  std::size_t lines = 0;
  for (std::size_t pos = 0; (pos = all_out.find("OUTPUT_LINE", pos)) != std::string::npos; pos += 11)
  {
    ++lines;
  }
  TEST_EQUAL(lines, 10000)
  const char* current_value = std::getenv("OPENMS_EXTERNAL_PROCESS_TEST");
  TEST_EQUAL(current_value != nullptr, inherited_present)
  TEST_EQUAL(current_value == nullptr ? std::string {} : std::string(current_value), original_value)
}
END_SECTION

START_SECTION([EXTRA] Windows environment overrides are case insensitive)
{
#ifdef OPENMS_WINDOWSPLATFORM
  const char* original = std::getenv("OPENMS_CASE_OVERRIDE_TEST");
  const std::string saved = original ? original : "";
  _putenv_s("OPENMS_CASE_OVERRIDE_TEST", "parent");
  std::string output;
  ExternalProcess process([&](const std::string& text) { output += text; }, [](const std::string&) {});
  const auto result = process.runWithResult(self_executable, {"--openms-process-environment", "OPENMS_CASE_OVERRIDE_TEST"}, "", false,
                                            ExternalProcess::IO_MODE::READ_WRITE, {{"openms_case_override_test", "child"}});
  _putenv_s("OPENMS_CASE_OVERRIDE_TEST", saved.c_str());
  TEST_EQUAL(result.state, ExternalProcess::RETURNSTATE::SUCCESS)
  TEST_EQUAL(output, "child")
#else
  NOT_TESTABLE
#endif
}
END_SECTION

START_SECTION([EXTRA] failed or crashed direct children cannot leave descendants writing outputs)
{
#ifndef OPENMS_WINDOWSPLATFORM
  TempDir tmp_dir;
  ExternalProcess process;
  for (const auto& ending : {"exit 37", "kill -KILL $$"})
  {
    const auto marker = to_path(tmp_dir.getPath()) / "survived.txt";
    const auto result
      = process.runWithResult("/bin/sh", {"-c", std::string("(sleep 0.3; echo SURVIVED > survived.txt) & ") + ending}, tmp_dir.getPath(), false);
    TEST_NOT_EQUAL(result.state, ExternalProcess::RETURNSTATE::SUCCESS)
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    TEST_FALSE(std::filesystem::exists(marker))
    std::filesystem::remove(marker);
  }
#else
  NOT_TESTABLE
#endif
}
END_SECTION

START_SECTION([EXTRA] cancellation terminates the process group and drains inherited pipes)
{
  TempDir tmp_dir;
  std::atomic_bool cancel {false};
  std::string all_out;
  ExternalProcess ep([&](const std::string& text) { all_out += text; }, [](const std::string&) {});
#ifdef OPENMS_WINDOWSPLATFORM
  const std::string shell = "cmd";
  const std::vector<std::string> long_args = {"/C", "echo READY & ping -n 30 127.0.0.1 >nul"};
#else
  const std::string shell = "/bin/sh";
  // READY comes from the descendant, ensuring it exists before cancellation.
  const std::vector<std::string> long_args = {"-c", "(echo READY; sleep 2; echo SURVIVED > descendant.txt) & wait"};
#endif
  const auto start = std::chrono::steady_clock::now();
  const auto result = ep.runWithResult(
    shell, long_args, tmp_dir.getPath(), false, ExternalProcess::IO_MODE::READ_WRITE, {},
    [&]() {
      if (all_out.find("READY") != std::string::npos || std::chrono::steady_clock::now() - start > std::chrono::seconds(3)) { cancel.store(true); }
    },
    &cancel);
  TEST_TRUE(all_out.find("READY") != std::string::npos)
  TEST_EQUAL(result.state, ExternalProcess::RETURNSTATE::CANCELLED)
  TEST_TRUE(std::chrono::steady_clock::now() - start < std::chrono::seconds(3))
#ifndef OPENMS_WINDOWSPLATFORM
  // Killing only the direct shell would leave its background descendant running.
  std::this_thread::sleep_for(std::chrono::milliseconds(2300));
  TEST_FALSE(std::filesystem::exists(std::filesystem::path(tmp_dir.getPath()) / "descendant.txt"))

  // Crash classification must also work when output is not captured.
  const auto crash = ep.runWithResult(shell, {"-c", "kill -KILL $$"}, "", false, ExternalProcess::IO_MODE::NO_IO);
  TEST_EQUAL(crash.state, ExternalProcess::RETURNSTATE::CRASH)
#endif
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#undef main

namespace
{
int dispatchProcessTest(int argc, char** argv)
{
  if (argc >= 2 && std::string(argv[1]) == "--openms-process-argv")
  {
    for (int index = 2; index < argc; ++index)
    {
      const std::string value(argv[index]);
      std::cout << value.size() << ':' << value;
    }
    return 0;
  }
  if (argc == 3 && std::string(argv[1]) == "--openms-process-environment")
  {
    const char* value = std::getenv(argv[2]);
    if (! value) return 1;
    std::cout << value;
    return 0;
  }
  if (argc == 3 && std::string(argv[1]) == "--openms-process-copy-env")
  {
#ifdef OPENMS_WINDOWSPLATFORM
    const wchar_t* destination = _wgetenv(L"OPENMS_UTF8_OUTPUT");
    if (! destination) return 1;
    std::filesystem::copy_file(to_path(argv[2]), std::filesystem::path(destination));
#else
    const char* destination = std::getenv("OPENMS_UTF8_OUTPUT");
    if (! destination) return 1;
    std::filesystem::copy_file(to_path(argv[2]), to_path(destination));
#endif
    return 0;
  }
  self_executable = pathString(std::filesystem::absolute(to_path(argv[0])));
  return externalProcessTestsMain(argc, argv);
}
} // namespace

#ifdef OPENMS_WINDOWSPLATFORM
// Use the CRT's wide entry point so Unicode checks validate what CreateProcessW
// delivered, without a second lossy conversion through the active code page.
int wmain(int argc, wchar_t** argv)
{
  std::vector<std::string> arguments;
  arguments.reserve(argc);
  for (int index = 0; index < argc; ++index)
    arguments.push_back(pathString(std::filesystem::path(argv[index])));
  std::vector<char*> pointers;
  for (auto& argument : arguments)
    pointers.push_back(argument.data());
  pointers.push_back(nullptr);
  return dispatchProcessTest(argc, pointers.data());
}
#else
int main(int argc, char** argv)
{ return dispatchProcessTest(argc, argv); }
#endif
