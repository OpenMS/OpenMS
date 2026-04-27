// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/InterProcessFileLock.h>

#include <chrono>
#include <thread>

#ifndef OPENMS_WINDOWSPLATFORM
#  include <signal.h>
#  include <sys/types.h>
#  include <sys/wait.h>
#  include <unistd.h>
#endif

///////////////////////////

using namespace OpenMS;

START_TEST(InterProcessFileLock, "$Id$")

// Acquire & release lock
START_SECTION((InterProcessFileLock(const String& target_file_path)))
{
  const String target = File::getTempDirectory() + "/ipfl_ctor_" + File::getUniqueName(false) + ".tmp";
  const String lock_file = target;

  {
    InterProcessFileLock lock(target);
    TEST_EQUAL(File::exists(lock_file), true)
  }

  TEST_EQUAL(File::remove(lock_file), true)
}
END_SECTION

// Re-acquire lock after release
START_SECTION((InterProcessFileLock(const String& target_file_path)))
{
  const String target = File::getTempDirectory() + "/ipfl_reacquire_" + File::getUniqueName(false) + ".tmp";
  const String lock_file = target;

  {
    InterProcessFileLock lock_a(target);
  }
  {
    InterProcessFileLock lock_b(target);
  }

  TEST_EQUAL(File::exists(lock_file), true)
  TEST_EQUAL(File::remove(lock_file), true)
}
END_SECTION

// validate constructor error handling for missing parent directory
START_SECTION((InterProcessFileLock(const String& target_file_path)))
{
  const String parent = File::getTempDirectory() + "/ipfl_missing_parent_" + File::getUniqueName(false);
  const String invalid_target = parent + "/child/target.tmp";
  TEST_EXCEPTION(Exception::FileNotWritable, InterProcessFileLock lock{invalid_target})
}
END_SECTION

// POSIX-specific test for lock contention and timeout behavior
START_SECTION((InterProcessFileLock(const String& target_file_path)))
{
#ifdef OPENMS_WINDOWSPLATFORM
  NOT_TESTABLE
#else
  const String target = File::getTempDirectory() + "/ipfl_contention_timeout_" + File::getUniqueName(false) + ".tmp";

  const pid_t child_pid = fork();
  TEST_NOT_EQUAL(child_pid, -1)

  if (child_pid == 0)
  {
    try
    {
      InterProcessFileLock child_lock(target);
      std::this_thread::sleep_for(std::chrono::seconds(20));
    }
    catch (const Exception::BaseException&)
    {
      _exit(1);
    }
    catch (const std::exception&)
    {
      _exit(1);
    }
    _exit(0);
  }

  std::this_thread::sleep_for(std::chrono::milliseconds(200));

  // Attempt to acquire the lock while it's held by the child process, expecting a timeout
  const auto start = std::chrono::steady_clock::now();
  bool timeout_thrown = false;
  try
  {
    InterProcessFileLock waiting_lock(target);
  }
  catch (const Exception::FailedAPICall&)
  {
    timeout_thrown = true;
  }

  const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::steady_clock::now() - start);

  TEST_TRUE(timeout_thrown)
  TEST_TRUE(elapsed >= std::chrono::milliseconds(2500))

  kill(child_pid, SIGTERM);
  int status = 0;
  waitpid(child_pid, &status, 0);

  TEST_EQUAL(File::remove(target), true)
#endif
}
END_SECTION

END_TEST
