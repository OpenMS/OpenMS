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
#include <OpenMS/CONCEPT/LogConfigHandler.h>
///////////////////////////

#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;

START_TEST(LogConfigHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((virtual ~LogConfigHandler()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Param parse(const StringList &setting)))
{
  StringList settings;
  settings.push_back("DEBUG add cout");
  settings.push_back("DEBUG add a.out");
  settings.push_back("INFO add a.out");
  settings.push_back("FATAL_ERROR add cerr");

  Param p = LogConfigHandler::getInstance()->parse(settings);

  // p should contain a list of the above set commands
  std::vector<std::string> parsedConfigs = p.getValue(LogConfigHandler::PARAM_NAME);

  TEST_EQUAL(parsedConfigs[0] , "DEBUG add cout FILE")
  TEST_EQUAL(parsedConfigs[1] , "DEBUG add a.out FILE")
  TEST_EQUAL(parsedConfigs[2] , "INFO add a.out FILE")
  TEST_EQUAL(parsedConfigs[3] , "FATAL_ERROR add cerr FILE")

  StringList settings2;
  settings2.push_back("DEBUG");

  TEST_EXCEPTION(Exception::ParseError, LogConfigHandler::getInstance()->parse(settings2));
}
END_SECTION

START_SECTION((void configure(const Param &param)))
{
  std::vector<std::string> settings = {"INFO add testing_info_warn_stream STRING",
                                      "WARNING add testing_info_warn_stream STRING",
                                      "ERROR add only_error_string_stream STRING",
                                      "INFO remove cout FILE",
                                      "WARNING remove cout",
                                      "ERROR remove cerr FILE"};

  Param p;
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");

  LogConfigHandler::getInstance()->configure(p);

  OPENMS_LOG_INFO << "1" << endl;
  OPENMS_LOG_INFO << "2" << endl;
  OPENMS_LOG_WARN << "3" << endl;
  OPENMS_LOG_ERROR << "4" << endl;

  settings.clear();
  settings.push_back("WARNING clear");
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");

  LogConfigHandler::getInstance()->configure(p);

  // this should go into nowhere
  OPENMS_LOG_WARN << "5" << endl;

  ostringstream& info_warn_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance()->getStream("testing_info_warn_stream"));
  String info_warn_stream_content(info_warn_stream.str());
  StringList info_warn_result;
  info_warn_stream_content.trim().split('\n', info_warn_result, true );

  TEST_EQUAL(info_warn_result.size() , 3)

  // check output with regex
  String pattern("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] ");
  boost::regex rx(pattern);

  int i = 1;
  for(StringList::iterator it = info_warn_result.begin() ; it != info_warn_result.end(); ++it)
  {
    rx.assign(pattern + i);
    TEST_TRUE(regex_search(*it, rx)) // stream may be wrapped in ANSI color codes; only search infix
    ++i;
  }
  ostringstream& error_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance()->getStream("only_error_string_stream"));
  String error_stream_content(error_stream.str());
  StringList error_result;
  error_stream_content.trim().split('\n', error_result, true );


  TEST_EQUAL(error_result.size(), 1)

  String pattern2("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] 4");
  rx.assign(pattern2);
  TEST_TRUE(regex_search(error_result[0], rx)) // stream may be wrapped in ANSI color codes; only search infix
}
END_SECTION

START_SECTION((ostream& getStream(const String &stream_name)))
{
  std::vector<std::string> settings;
  settings.push_back("INFO add testing_getStream STRING");

  Param p;
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");

  LogConfigHandler::getInstance()->configure(p);

  OPENMS_LOG_INFO << "getStream 1" << endl;

  ostringstream& info_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance()->getStream("testing_getStream"));
  String info_content(info_stream.str());

  StringList info_result;
  info_content.trim().split('\n', info_result, true );

  TEST_EQUAL(info_result.size() , 1)

  // check if everything landed in the stream we wanted
  String pattern("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] getStream 1");
  boost::regex rx(pattern);
  TEST_EQUAL(regex_match(info_result[0], rx), true)
}
END_SECTION

LogConfigHandler* nullPointer = nullptr;
START_SECTION((static LogConfigHandler* getInstance()))
{
  TEST_NOT_EQUAL(LogConfigHandler::getInstance(), nullPointer)
}
END_SECTION

START_SECTION((void setLogLevel(const String &log_level) - restoring streams))
{
  // Test that setLogLevel can restore streams when lowering the log level
  
  // Setup: Create a string stream for INFO level
  std::vector<std::string> settings;
  settings.push_back("INFO add test_setloglevel_stream STRING");
  
  Param p;
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");
  LogConfigHandler::getInstance()->configure(p);
  
  // Write a message at INFO level - should appear
  OPENMS_LOG_INFO << "message1" << endl;
  
  // Set log level to ERROR (should remove INFO streams)
  LogConfigHandler::getInstance()->setLogLevel("ERROR");
  
  // Write a message at INFO level - should NOT appear
  OPENMS_LOG_INFO << "message2" << endl;
  
  // Lower log level back to INFO (should restore INFO streams)
  LogConfigHandler::getInstance()->setLogLevel("INFO");
  
  // Write a message at INFO level - should appear again
  OPENMS_LOG_INFO << "message3" << endl;
  
  // Check the stream content
  ostringstream& test_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance()->getStream("test_setloglevel_stream"));
  String content(test_stream.str());
  StringList result;
  content.trim().split('\n', result, true);
  
  // Should have message1 and message3, but not message2
  TEST_EQUAL(result.size(), 2)
  TEST_TRUE(result[0].hasSubstring("message1"))
  TEST_TRUE(result[1].hasSubstring("message3"))
}
END_SECTION

START_SECTION((removeAllStreams flushes buffers))
{
  // Test that removeAllStreams flushes buffers before clearing
  
  // Setup: Create a string stream for WARNING level
  std::vector<std::string> settings;
  settings.push_back("WARNING add test_flush_stream STRING");
  
  Param p;
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");
  LogConfigHandler::getInstance()->configure(p);
  
  // Write without endl (no flush yet)
  OPENMS_LOG_WARN << "unflushed_message";
  
  // Set log level to ERROR (which calls removeAllStreams on WARNING)
  // This should flush the buffer before removing streams
  LogConfigHandler::getInstance()->setLogLevel("ERROR");
  
  // Check the stream content - the unflushed message should be there
  ostringstream& test_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance()->getStream("test_flush_stream"));
  String content(test_stream.str());
  
  // The message should be in the stream even though it wasn't flushed with endl
  TEST_TRUE(content.hasSubstring("unflushed_message"))
}
END_SECTION

START_SECTION((void setLogLevel(const String &log_level) - NONE level))
{
  // Test that setLogLevel("NONE") disables all logging and can be restored
  
  // Setup: Create string streams for multiple levels
  std::vector<std::string> settings;
  settings.push_back("INFO add test_none_info_stream STRING");
  settings.push_back("ERROR add test_none_error_stream STRING");
  
  Param p;
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");
  LogConfigHandler::getInstance()->configure(p);
  
  // Write messages - should appear
  OPENMS_LOG_INFO << "before_none_info" << endl;
  OPENMS_LOG_ERROR << "before_none_error" << endl;
  
  // Set log level to NONE (should remove all streams)
  LogConfigHandler::getInstance()->setLogLevel("NONE");
  
  // Write messages - should NOT appear
  OPENMS_LOG_INFO << "during_none_info" << endl;
  OPENMS_LOG_ERROR << "during_none_error" << endl;
  
  // Restore log level to INFO (should restore INFO and higher streams)
  LogConfigHandler::getInstance()->setLogLevel("INFO");
  
  // Write messages - should appear again
  OPENMS_LOG_INFO << "after_none_info" << endl;
  OPENMS_LOG_ERROR << "after_none_error" << endl;
  
  // Check INFO stream
  ostringstream& info_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance()->getStream("test_none_info_stream"));
  String info_content(info_stream.str());
  TEST_TRUE(info_content.hasSubstring("before_none_info"))
  TEST_FALSE(info_content.hasSubstring("during_none_info"))
  TEST_TRUE(info_content.hasSubstring("after_none_info"))
  
  // Check ERROR stream
  ostringstream& error_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance()->getStream("test_none_error_stream"));
  String error_content(error_stream.str());
  TEST_TRUE(error_content.hasSubstring("before_none_error"))
  TEST_FALSE(error_content.hasSubstring("during_none_error"))
  TEST_TRUE(error_content.hasSubstring("after_none_error"))
}
END_SECTION

START_SECTION((void setLogLevel(const String &log_level) - invalid level))
{
  // Test that setLogLevel throws an exception for invalid log levels
  TEST_EXCEPTION(Exception::IllegalArgument, LogConfigHandler::getInstance()->setLogLevel("INVALID_LEVEL"))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
