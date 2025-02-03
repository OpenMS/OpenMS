// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/StreamHandler.h>
///////////////////////////

#include <sstream>

using namespace OpenMS;
using namespace std;

START_TEST(StreamHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

StreamHandler* ptr = nullptr;
StreamHandler* nullPointer = nullptr;
START_SECTION(StreamHandler())
{
	ptr = new StreamHandler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~StreamHandler())
{
	delete ptr;
}
END_SECTION

// main instance for the test
StreamHandler handler;

START_SECTION((Int registerStream(StreamType const type, const String &stream_name)))
{
  String filename;
  NEW_TMP_FILE(filename)

  handler.registerStream(StreamHandler::FILE, filename);

  handler.getStream(StreamHandler::FILE, filename) << "This is a test!" << endl;

  ostream & s = handler.getStream(StreamHandler::FILE, filename);

  s << "And another test!" << endl;

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("StreamHandler_test.txt"))

  // if you try to register a stream with the same name, but a different type
  // an Exception should be thrown
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, handler.registerStream(StreamHandler::STRING, filename), "This stream was already registered with a different type.")
}
END_SECTION

START_SECTION((void unregisterStream(StreamType const type, const String &stream_name)))
{
  String filename;
  NEW_TMP_FILE(filename)

  // this one was registered twice
  handler.registerStream(StreamHandler::FILE, filename);
  handler.registerStream(StreamHandler::FILE, filename);

  // one unregister .. it should still be available
  handler.unregisterStream(StreamHandler::FILE, filename);

  handler.getStream(StreamHandler::FILE, filename);

  // now it should be gone
  handler.unregisterStream(StreamHandler::FILE, filename);

  TEST_EXCEPTION(Exception::ElementNotFound, handler.unregisterStream(StreamHandler::FILE, filename))
}
END_SECTION

START_SECTION((ostream& getStream(StreamType const type, const String &stream_name)))
{
  String file2;
  NEW_TMP_FILE(file2);

  handler.registerStream(StreamHandler::FILE, file2);
  handler.getStream(StreamHandler::FILE, file2) << "This is a test!" << endl;

  ostream & file_stream = handler.getStream(StreamHandler::FILE, file2);
  file_stream << "And another test!" << endl;

  TEST_FILE_EQUAL(file2.c_str(), OPENMS_GET_TEST_DATA_PATH("StreamHandler_test.txt"))

  // now we test this with stringstreams
  handler.registerStream(StreamHandler::STRING, "getStream_testing_stream");
  handler.getStream(StreamHandler::STRING, "getStream_testing_stream") << "This is a test!" << endl;

  ostream & string_stream = handler.getStream(StreamHandler::STRING, "getStream_testing_stream");
  string_stream << "And another test!" << endl;

  ostringstream & ostr = static_cast<ostringstream&>(handler.getStream(StreamHandler::STRING, "getStream_testing_stream"));
  String output(ostr.str());
  StringList results;
  output.trim().split('\n',results);

  TEST_EQUAL(results.size(), 2)
  TEST_EQUAL(results[0], "This is a test!")
  TEST_EQUAL(results[1], "And another test!")
}
END_SECTION

START_SECTION((bool hasStream(const StreamType type, const String &stream_name)))
{
  handler.registerStream(StreamHandler::STRING, "this_is_a_test_stream");

  TEST_EQUAL(handler.hasStream(StreamHandler::STRING, "this_is_a_test_stream"), true)
  TEST_EQUAL(handler.hasStream(StreamHandler::FILE, "this_is_a_test_stream"), false)
  TEST_EQUAL(handler.hasStream(StreamHandler::STRING, "this_is_not_the_same_stream"), false)
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
