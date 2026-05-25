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
  // Create a unique lock file path in temp directory for isolated test execution
  const String target = File::getTempDirectory() + "/ipfl_ctor_" + File::getUniqueName(false) + ".tmp";
  const String lock_file = target;

  {
    // Constructor should acquire lock immediately and create lock file
    InterProcessFileLock lock(target);
    TEST_TRUE(lock.isLocked())
    TEST_TRUE(File::exists(lock_file))
  }

  // After RAII object destruction, file should still exist and be removable
  TEST_TRUE(File::remove(lock_file))
}
END_SECTION

// Re-acquire lock after release
START_SECTION((InterProcessFileLock(const String& target_file_path)))
{
  // Use a second unique path to avoid cross-test interference
  const String target = File::getTempDirectory() + "/ipfl_reacquire_" + File::getUniqueName(false) + ".tmp";
  const String lock_file = target;

  {
    // First lock acquisition should succeed
    InterProcessFileLock lock_a(target);
    TEST_TRUE(lock_a.isLocked())
  }
  {
    // Re-acquisition on same path should succeed after previous lock scope ends
    InterProcessFileLock lock_b(target);
    TEST_TRUE(lock_b.isLocked())
  }

  // Ensure file is still present and can be cleaned up
  TEST_TRUE(File::exists(lock_file))
  TEST_TRUE(File::remove(lock_file))
}
END_SECTION

// Validate lock() failure then recovery with a valid target
START_SECTION((bool timedLock(const String& target_file_path, UInt timeout_ms = 100) noexcept))
{
  // Build invalid path (missing parent directory tree) and one valid path
  const String parent = File::getTempDirectory() + "/ipfl_missing_parent_" + File::getUniqueName(false);
  const String invalid_target = parent + "/child/target.tmp";
  const String valid_target = File::getTempDirectory() + "/ipfl_lock_valid_" + File::getUniqueName(false) + ".tmp";

  {
    // API is non-throwing: failure is reported via return value and unlocked state
    InterProcessFileLock lock;
    TEST_FALSE(lock.timedLock(invalid_target, 1000))
    TEST_FALSE(lock.isLocked())

    // Subsequent timedLock() call with valid target should recover and lock
    TEST_TRUE(lock.timedLock(valid_target, 1000))
    TEST_TRUE(lock.isLocked())
    TEST_TRUE(File::exists(valid_target))
  }

  // Cleanup
  TEST_TRUE(File::remove(valid_target))
}
END_SECTION

START_SECTION((bool tryLock(const String& target_file_path) noexcept))
{
  const String parent = File::getTempDirectory() + "/ipfl_trylock_missing_parent_" + File::getUniqueName(false);
  const String invalid_target = parent + "/child/target.tmp";

  InterProcessFileLock lock;
  TEST_FALSE(lock.tryLock(invalid_target))
  TEST_FALSE(lock.isLocked())
}
END_SECTION

// Validate timeout_ms == 0 normalization (implementation enforces minimum of 1ms)
START_SECTION((InterProcessFileLock(const String& target_file_path, UInt timeout_ms)))
{
  // Lock with zero timeout must still succeed on uncontended target
  const String target = File::getTempDirectory() + "/ipfl_timeout_zero_" + File::getUniqueName(false) + ".tmp";

  {
    InterProcessFileLock lock(target, 0);
    TEST_TRUE(lock.isLocked())
    TEST_TRUE(File::exists(target))
  }

  // Cleanup lock file
  TEST_TRUE(File::remove(target))
}
END_SECTION

// validate constructor error handling for missing parent directory (non-throwing)
START_SECTION((InterProcessFileLock(const String& target_file_path)))
{
  // Build a deliberately invalid path (missing parent directory tree)
  const String parent = File::getTempDirectory() + "/ipfl_missing_parent_" + File::getUniqueName(false);
  const String invalid_target = parent + "/child/target.tmp";

  // API is non-throwing: failure is reported via unlocked state
  InterProcessFileLock lock{invalid_target};
  TEST_FALSE(lock.isLocked())
}
END_SECTION

// POSIX-specific test for lock contention and timeout behavior
START_SECTION((InterProcessFileLock(const String& target_file_path, UInt timeout_ms)))
{
#ifdef OPENMS_WINDOWSPLATFORM
  NOT_TESTABLE
#else
  // Shared target lock file; child acquires lock first, parent tries to acquire with timeout
  const String target = File::getTempDirectory() + "/ipfl_contention_timeout_" + File::getUniqueName(false) + ".tmp";
  const UInt timeout_ms = 1000;

  // One-byte pipe acts as a readiness barrier (child -> parent)
  int ready_pipe[2];
  TEST_EQUAL(pipe(ready_pipe), 0)
  const int byte_amount = 1;
  const int read_end = 0;
  const int write_end = 1;

  // Fork a child process that holds the lock
  const pid_t child_pid = fork();
  TEST_NOT_EQUAL(child_pid, -1)

  if (child_pid == 0)
  {
    // Child only writes readiness; close read-end
    close(ready_pipe[0]);

    // Watchdog for the child process in case parent cleanup/termination fails
    const unsigned int child_termination_watchdog_seconds = static_cast<unsigned int>(timeout_ms / 1000 + 10);
    // schedule SIGALARM for the child after x seconds
    alarm(child_termination_watchdog_seconds);
    try
    {
      InterProcessFileLock child_lock(target, timeout_ms);

      // Notify parent that lock is held
      const char ok = '1';
      write(ready_pipe[write_end], &ok, byte_amount);

      // Wait for parent termination signal; watchdog above bounds runtime
      pause();
    }
    catch (const Exception::BaseException&)
    {
      const char fail = '0';
      write(ready_pipe[write_end], &fail, byte_amount);
      pause();
    }
    catch (const std::exception&)
    {
      const char fail = '0';
      write(ready_pipe[write_end], &fail, byte_amount);
      write(ready_pipe[write_end], &fail, byte_amount);
      pause();
    }
  }
  else
  {
    // Parent only reads readiness; close write-end
    close(ready_pipe[write_end]);

    // Wait until child reports lock acquisition success/failure
    char ready = 0;
    ssize_t bytes = read(ready_pipe[read_end], &ready, byte_amount);
    close(ready_pipe[read_end]);
    TEST_EQUAL(bytes, byte_amount)
    TEST_EQUAL(ready, '1')

    // Attempt to acquire the lock while it's held by the child process, expecting timeout/failure
    const auto start = std::chrono::steady_clock::now();
    InterProcessFileLock waiting_lock(target, timeout_ms);

    const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - start);

    TEST_FALSE(waiting_lock.isLocked())

    // The timeout should be at least the specified duration
    TEST_TRUE(elapsed >= std::chrono::seconds(1))

    // Ensure child process is terminated
    kill(child_pid, SIGTERM);
    int status = 0;
    waitpid(child_pid, &status, 0);

    // Cleanup lock file generated by test
    TEST_EQUAL(File::remove(target), true)
  }
#endif
}
END_SECTION

END_TEST
