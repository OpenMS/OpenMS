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
#include <string>
#include <thread>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Exception;
using namespace std;

START_TEST(GlobalExceptionHandler, "$Id$")

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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



