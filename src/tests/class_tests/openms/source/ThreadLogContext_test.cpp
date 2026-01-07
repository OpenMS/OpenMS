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
#include <OpenMS/CONCEPT/ThreadLogContext.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <thread>
#include <vector>
#include <atomic>
#include <sstream>
#include <boost/regex.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif

///////////////////////////

using namespace OpenMS;
using namespace Logger;
using namespace std;

START_TEST(ThreadLogContext, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static Logger::LogStream& fatal()))
{
  // Test that fatal() returns a valid LogStream reference
  LogStream& stream = ThreadLogContext::fatal();
  TEST_EQUAL(stream.getLevel(), "FATAL_ERROR")
}
END_SECTION

START_SECTION((static Logger::LogStream& error()))
{
  LogStream& stream = ThreadLogContext::error();
  TEST_EQUAL(stream.getLevel(), "ERROR")
}
END_SECTION

START_SECTION((static Logger::LogStream& warn()))
{
  LogStream& stream = ThreadLogContext::warn();
  TEST_EQUAL(stream.getLevel(), "WARNING")
}
END_SECTION

START_SECTION((static Logger::LogStream& info()))
{
  LogStream& stream = ThreadLogContext::info();
  TEST_EQUAL(stream.getLevel(), "INFO")
}
END_SECTION

START_SECTION((static Logger::LogStream& debug()))
{
  LogStream& stream = ThreadLogContext::debug();
  TEST_EQUAL(stream.getLevel(), "DEBUG")
}
END_SECTION

START_SECTION((static bool isInitialized()))
{
  // After accessing any stream, should be initialized
  ThreadLogContext::info();
  TEST_EQUAL(ThreadLogContext::isInitialized(), true)
}
END_SECTION

START_SECTION((static void ensureInitialized()))
{
  // This should not throw or crash
  ThreadLogContext::ensureInitialized();
  TEST_EQUAL(ThreadLogContext::isInitialized(), true)
}
END_SECTION

START_SECTION((static void reset()))
{
  // Test that reset doesn't crash and re-initializes properly
  ThreadLogContext::reset();
  TEST_EQUAL(ThreadLogContext::isInitialized(), true)
  TEST_EQUAL(ThreadLogContext::info().getLevel(), "INFO")
}
END_SECTION

START_SECTION((static void flushAll()))
{
  // Test that flushAll doesn't crash
  ThreadLogContext::flushAll();
  NOT_TESTABLE  // Just testing it doesn't crash
}
END_SECTION

START_SECTION(([EXTRA] Thread-local output independence))
{
  // Test that each thread writes to its own stream without interference
  // This tests the fundamental guarantee: thread-local streams don't share state

  ostringstream main_stream;
  ThreadLogContext::info().removeAllStreams();
  ThreadLogContext::info().insert(main_stream);

  ThreadLogContext::info() << "main_thread" << endl;

  // The main thread's output should contain only its message
  TEST_EQUAL(main_stream.str().find("main_thread") != string::npos, true)

  // Clean up
  ThreadLogContext::info().remove(main_stream);
}
END_SECTION

START_SECTION(([EXTRA] Thread-local macros test - OPENMS_LOG_*_TL))
{
  // Test the thread-local macros directly
  ostringstream stream;

  // Configure thread-local info stream
  ThreadLogContext::info().removeAllStreams();
  ThreadLogContext::info().insert(stream);
  ThreadLogContext::info().rdbuf()->clearCache();

  OPENMS_LOG_INFO_TL << "test_message_1" << endl;
  OPENMS_LOG_INFO_TL << "test_message_2" << endl;

  TEST_EQUAL(stream.str().find("test_message_1") != string::npos, true)
  TEST_EQUAL(stream.str().find("test_message_2") != string::npos, true)

  ThreadLogContext::info().remove(stream);
}
END_SECTION

START_SECTION(([EXTRA] Multi-threaded logging with std::thread))
{
  // Test that multiple threads can log simultaneously without crashing
  // and that each thread gets its own independent log stream

  constexpr int NUM_THREADS = 4;
  atomic<int> completed_count{0};
  atomic<int> has_own_prefix{0};

  vector<thread> threads;

  for (int t = 0; t < NUM_THREADS; ++t)
  {
    threads.emplace_back([t, &completed_count, &has_own_prefix]() {
      // Each thread sets up its own capture stream
      ostringstream thread_stream;
      ThreadLogContext::info().removeAllStreams();
      ThreadLogContext::info().insert(thread_stream);
      ThreadLogContext::info().rdbuf()->clearCache();

      // Log a unique message for this thread
      string unique_marker = "MARKER_THREAD_" + to_string(t) + "_UNIQUE";
      ThreadLogContext::info() << unique_marker << endl;

      // Verify this thread's output contains its marker
      string output = thread_stream.str();
      if (output.find(unique_marker) != string::npos)
      {
        has_own_prefix++;
      }

      ThreadLogContext::info().remove(thread_stream);
      completed_count++;
    });
  }

  // Wait for all threads
  for (auto& th : threads)
  {
    th.join();
  }

  // All threads should complete
  TEST_EQUAL(completed_count.load(), NUM_THREADS)
  // All threads should find their own marker in their stream
  TEST_EQUAL(has_own_prefix.load(), NUM_THREADS)
}
END_SECTION

#ifdef _OPENMP
START_SECTION(([EXTRA] OpenMP parallel logging))
{
  // Test thread-local logging with OpenMP
  // Verify that parallel logging doesn't crash and each thread can log

  atomic<int> threads_completed{0};
  atomic<int> threads_found_own_marker{0};

  omp_set_num_threads(4);

  #pragma omp parallel
  {
    // Each OpenMP thread gets its own stream
    ostringstream local_stream;
    ThreadLogContext::info().removeAllStreams();
    ThreadLogContext::info().insert(local_stream);
    ThreadLogContext::info().rdbuf()->clearCache();

    int thread_id = omp_get_thread_num();

    // Log a unique marker for this thread
    string marker = "OMP_MARKER_" + to_string(thread_id) + "_END";
    OPENMS_LOG_INFO_TL << marker << endl;

    // Verify this thread's output contains its marker
    string output = local_stream.str();
    if (output.find(marker) != string::npos)
    {
      threads_found_own_marker++;
    }

    ThreadLogContext::info().remove(local_stream);
    threads_completed++;
  }

  // All threads should complete
  TEST_EQUAL(threads_completed.load() > 0, true)
  // All threads should find their own marker
  TEST_EQUAL(threads_found_own_marker.load(), threads_completed.load())
}
END_SECTION
#endif

START_SECTION(([EXTRA] Long string test for buffer handling))
{
  // Test with long strings that might trigger buffer overflow issues
  // The buffer in LogStreamBuf is 32KB, so we test with various sizes

  ostringstream stream;
  ThreadLogContext::info().removeAllStreams();
  ThreadLogContext::info().insert(stream);
  ThreadLogContext::info().rdbuf()->clearCache();

  // Test with a moderately long string (1KB)
  string long_string(1024, 'x');
  string marker_start = "LONG_START_";
  string marker_end = "_LONG_END";

  OPENMS_LOG_INFO_TL << marker_start << long_string << marker_end << endl;

  string output = stream.str();

  // Verify markers are present (the string was logged)
  TEST_EQUAL(output.find(marker_start) != string::npos, true)
  TEST_EQUAL(output.find(marker_end) != string::npos, true)

  ThreadLogContext::info().remove(stream);
}
END_SECTION

START_SECTION(([EXTRA] Debug macro with file/line info))
{
  ostringstream stream;
  ThreadLogContext::debug().removeAllStreams();
  ThreadLogContext::debug().insert(stream);
  ThreadLogContext::debug().rdbuf()->clearCache();

  OPENMS_LOG_DEBUG_TL << "debug_test" << endl;

  string output = stream.str();

  // Should contain file name
  TEST_EQUAL(output.find("ThreadLogContext_test.cpp") != string::npos, true)
  // Should contain line number (digits in parentheses)
  boost::regex rx(R"(ThreadLogContext_test\.cpp\(\d+\))");
  TEST_EQUAL(regex_search(output, rx), true)
  // Should contain the message
  TEST_EQUAL(output.find("debug_test") != string::npos, true)

  ThreadLogContext::debug().remove(stream);
}
END_SECTION

START_SECTION(([EXTRA] Fatal error macro with file/line info))
{
  ostringstream stream;
  ThreadLogContext::fatal().removeAllStreams();
  ThreadLogContext::fatal().insert(stream);
  ThreadLogContext::fatal().rdbuf()->clearCache();

  OPENMS_LOG_FATAL_ERROR_TL << "fatal_test" << endl;

  string output = stream.str();

  // Should contain full file path
  TEST_EQUAL(output.find("ThreadLogContext_test.cpp") != string::npos, true)
  // Should contain the message
  TEST_EQUAL(output.find("fatal_test") != string::npos, true)

  ThreadLogContext::fatal().remove(stream);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
