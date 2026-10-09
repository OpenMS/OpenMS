// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>
///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>

#include <atomic>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// The terminate() section starts this test program as a child process (by /proc/self/exe; Linux only).
#if defined(__linux__)
#define TERMINATE_CHILD_TESTS 1
#include <cerrno>
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#else
#define TERMINATE_CHILD_TESTS 0
#endif

using namespace OpenMS;
using namespace OpenMS::Exception;
using namespace std;

START_TEST(GlobalExceptionHandler, "$Id$")

#if TERMINATE_CHILD_TESTS
// child process of the terminate() section below: an exception constructed on a worker thread and rethrown, not
// caught, on this one (as ProSE's parallel loops rethrow the first exception of a worker thread)
if (std::getenv("OPENMS_GEH_TEST_TERMINATE_CHILD") != nullptr)
{
  std::exception_ptr transferred;
  std::thread worker([&transferred]()
  {
    try
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, "worker()", "the transferred exception", "42");
    }
    catch (...)
    {
      transferred = std::current_exception();
    }
  });
  worker.join();
  // not caught (START_TEST catches exceptions, a noexcept function does not let them pass): terminate()
  [&transferred]() noexcept { std::rethrow_exception(transferred); }();
}
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static GlobalExceptionHandler& getInstance()))
{
  // the singleton; its constructor installs the terminate and new handlers
  GlobalExceptionHandler& instance = GlobalExceptionHandler::getInstance();
  TEST_EQUAL(&instance == &GlobalExceptionHandler::getInstance(), true)
}
END_SECTION

START_SECTION((static void set(const std::string& file, int line, const std::string& function, const std::string& name, const std::string& message)))
{
  GlobalExceptionHandler::set("file", 1, "function", "name", "message");
  GlobalExceptionHandler::setName("name");
  GlobalExceptionHandler::setMessage("message");
  GlobalExceptionHandler::setLine(2);
  GlobalExceptionHandler::setFile("file");
  GlobalExceptionHandler::setFunction("function");
  NOT_TESTABLE // the state is read by terminate() only
}
END_SECTION

START_SECTION(([EXTRA] exceptions constructed on several threads at once record their state per thread))
{
  // Every exception constructor records file, line, function, name and message in the handler. Threads
  // that throw at the same time (e.g. the chunks of a parallel file read that fail to parse) used to
  // write the same strings concurrently (a data race ThreadSanitizer reports as a heap-use-after-free);
  // the state is thread-local now. This runs the pattern; the race itself needs a sanitizer to show.
  std::vector<std::thread> threads;
  std::atomic<int> caught{0};
  for (int t = 0; t < 4; ++t)
  {
    threads.emplace_back([t, &caught]()
    {
      for (int i = 0; i < 1000; ++i)
      {
        try
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "thread " + std::to_string(t) + " iteration " + std::to_string(i) + std::string(64, 'x'));
        }
        catch (const Exception::BaseException&)
        {
          ++caught;
        }
      }
    });
  }
  for (auto& thread : threads) thread.join();
  TEST_EQUAL(caught.load(), 4000)
}
END_SECTION


START_SECTION((static void terminate()))
{
#if TERMINATE_CHILD_TESTS
  // The state is recorded per thread: for an exception that another thread constructed (and std::exception_ptr
  // carried here), the last entry of this thread is about another exception or none. terminate() reports the
  // exception that was not caught.
  if (::access("/proc/self/exe", X_OK) != 0)
  {
    STATUS("SKIPPED: /proc/self/exe is not available, so this test program cannot start itself as the child process")
  }
  else
  {
    std::string out;
    NEW_TMP_FILE(out)
    std::vector<std::string> env_strings;
    for (char** e = environ; *e != nullptr; ++e) env_strings.emplace_back(*e);
    env_strings.emplace_back("OPENMS_GEH_TEST_TERMINATE_CHILD=1");
    std::vector<char*> envp;
    for (std::string& e : env_strings) envp.push_back(&e[0]);
    envp.push_back(nullptr);
    char arg0[] = "GlobalExceptionHandler_test";
    char* child_argv[] = {arg0, nullptr};
    posix_spawn_file_actions_t actions;
    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, out.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
    pid_t pid = 0;
    const int spawn_rc = ::posix_spawn(&pid, "/proc/self/exe", &actions, nullptr, child_argv, envp.data());
    posix_spawn_file_actions_destroy(&actions);
    TEST_EQUAL(spawn_rc, 0)
    int status = 0;
    pid_t waited = -1;
    if (spawn_rc == 0)
    {
      do { waited = ::waitpid(pid, &status, 0); } while (waited < 0 && errno == EINTR);
    }
    TEST_EQUAL(waited == pid && WIFSIGNALED(status), true) // abort()
    std::ifstream in(out);
    std::stringstream report;
    report << in.rdbuf();
    TEST_EQUAL(report.str().find("FATAL: uncaught exception!") != std::string::npos, true)
    TEST_EQUAL(report.str().find("the transferred exception") != std::string::npos, true)
    TEST_EQUAL(report.str().find("InvalidValue") != std::string::npos, true)
  }
#else
  STATUS("SKIPPED: the child process is started by /proc/self/exe (Linux only)")
#endif
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



