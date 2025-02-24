// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Stephan Aiche $
// $Authors: Chris Bielow, Stephan Aiche, Andreas Bertsch $
// --------------------------------------------------------------------------


/**

  Most of the tests, generously provided by the BALL people, taken from version 1.2

*/

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <fstream>
#include <boost/regex.hpp>

// OpenMP support
#ifdef _OPENMP
	#include <omp.h>
#endif


///////////////////////////

using namespace OpenMS;
using namespace Logger;
using namespace std;

class TestTarget
  :  public LogStreamNotifier
{
  public:
  void logNotify() override
  {
    notified = true;
    return;
  }
  bool notified;
};

START_TEST(LogStream, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] OpenMP - test))
{
  // just see if this crashes with OpenMP
  ostringstream stream_by_logger;
  OpenMS_Log_debug.insert(stream_by_logger);
  OpenMS_Log_debug.remove(cout);
  OpenMS_Log_info.insert(stream_by_logger);
  OpenMS_Log_info.remove(cout);

  {
    // create a long string that is of similar length as the buffer length to
    // ensure buffering and flushing works correctly LogStream.cpp even in a
    // multi-threaded environment.
    std::string long_str;
    for (int k = 0; k < 32768/2; k++) if (char(k) != 0) long_str += char(k);

    #ifdef _OPENMP
    omp_set_num_threads(8);
    #pragma omp parallel for
    #endif
    for (int i=0;i<10000;++i)
    {
      OPENMS_LOG_DEBUG << long_str << "1\n";
      OPENMS_LOG_DEBUG << "2" << endl;
      OPENMS_LOG_INFO << "1\n";
      OPENMS_LOG_INFO << "2" << endl;
    }
  }

  // remove logger after testing
  OpenMS_Log_debug.remove(stream_by_logger);
  OpenMS_Log_info.remove(stream_by_logger);

  NOT_TESTABLE;
}
END_SECTION

LogStream* nullPointer = nullptr;

START_SECTION(LogStream(LogStreamBuf *buf=0, bool delete_buf=true, std::ostream* stream))
{
  LogStream* l1 = new LogStream((LogStreamBuf*)nullptr);
  TEST_NOT_EQUAL(l1, nullPointer)
  delete l1;

  LogStreamBuf* lb2(new LogStreamBuf());
  LogStream* l2 = new LogStream(lb2);
  TEST_NOT_EQUAL(l2, nullPointer)
  delete l2;
}
END_SECTION

START_SECTION((virtual ~LogStream()))
{
	ostringstream stream_by_logger;
  {
		LogStream* l1 = new LogStream(new LogStreamBuf());
		l1->insert(stream_by_logger);
		*l1 << "flushtest" << endl;
		TEST_EQUAL(stream_by_logger.str(),"flushtest\n")
		*l1 << "unfinishedline...";
		TEST_EQUAL(stream_by_logger.str(),"flushtest\n")
		delete l1;
		// testing if loggers' d'tor will distribute the unfinished line to its children...
	}
	TEST_EQUAL(stream_by_logger.str(),"flushtest\nunfinishedline...\n")

}
END_SECTION


START_SECTION((LogStreamBuf* operator->()))
{
  LogStream l1(new LogStreamBuf());
  l1->sync(); // if it doesn't crash we're happy
  NOT_TESTABLE
}
END_SECTION

START_SECTION((LogStreamBuf* rdbuf()))
{
  LogStream l1(new LogStreamBuf());
  // small workaround since TEST_NOT_EQUAL(l1.rdbuf, 0) would expand to
  // cout << ls.rdbuf()
  // which kills the cout buffer
  TEST_NOT_EQUAL((l1.rdbuf()==nullptr), true)
}
END_SECTION

START_SECTION((void setLevel(std::string level)))
{
  LogStream l1(new LogStreamBuf());
  l1.setLevel("INFORMATION");
  TEST_EQUAL(l1.getLevel(), "INFORMATION")
}
END_SECTION

START_SECTION((std::string getLevel()))
{
  LogStream l1(new LogStreamBuf());
  TEST_EQUAL(l1.getLevel(), LogStreamBuf::UNKNOWN_LOG_LEVEL)
  l1.setLevel("FATAL_ERROR");
  TEST_EQUAL(l1.getLevel(), "FATAL_ERROR")
}
END_SECTION

START_SECTION((void insert(std::ostream &s)))
{
  String filename;
  NEW_TMP_FILE(filename)
  LogStream l1(new LogStreamBuf());
  ofstream s(filename.c_str(), std::ios::out);
  l1.insert(s);

  l1 << "1\n";
  l1 << "2" << endl;

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION((void remove(std::ostream &s)))
{
  LogStream l1(new LogStreamBuf());
  ostringstream s;
  l1 << "BLA"<<endl;
  l1.insert(s);
  l1 << "to_stream"<<endl;
  l1.remove(s);
  // make sure we can remove it twice without harm
  l1.remove(s);
	l1 << "BLA2"<<endl;
  TEST_EQUAL(s.str(),"to_stream\n");
}
END_SECTION

START_SECTION((void insertNotification(std::ostream &s, LogStreamNotifier &target)))
{
  LogStream l1(new LogStreamBuf());
  TestTarget target;
  ofstream os;
  target.registerAt(l1);
  target.notified = false;
  TEST_EQUAL(target.notified, false)
  l1 << "test" << std::endl;
  TEST_EQUAL(target.notified, true)
}
END_SECTION

START_SECTION(([EXTRA]removeNotification))
{
  LogStream l1(new LogStreamBuf());
  TestTarget target;
  ofstream os;
  target.registerAt(l1);
  target.unregister();
  target.notified = false;
  TEST_EQUAL(target.notified, false)
  l1 << "test" << endl;
  TEST_EQUAL(target.notified, false)
  // make sure we can remove it twice
  target.unregister();
  l1 << "test" << endl;
  TEST_EQUAL(target.notified, false)
}
END_SECTION

START_SECTION((void setPrefix(const std::string &prefix)))
{
	LogStream l1(new LogStreamBuf());
	ostringstream stream_by_logger;
	l1.insert(stream_by_logger);
	l1.setLevel("DEVELOPMENT");
	l1.setPrefix("%y"); //message type ("Error", "Warning", "Information", "-")
	l1 << "  2." << endl;
	l1.setPrefix("%T"); //time (HH:MM:SS)
	l1 << "  3." << endl;
	l1.setPrefix( "%t"); //time in short format (HH:MM)
	l1 << "  4." << endl;
	l1.setPrefix("%D"); //date (YYYY/MM/DD)
	l1 << "  5." << endl;
	l1.setPrefix("%d"); // date in short format (MM/DD)
	l1 << "  6." << endl;
	l1.setPrefix("%S"); //time and date (YYYY/MM/DD, HH:MM:SS)
	l1 << "  7." << endl;
	l1.setPrefix("%s"); //time and date in short format (MM/DD, HH:MM)
	l1 << "  8." << endl;
	l1.setPrefix("%%"); //percent sign (escape sequence)
	l1 << "  9." << endl;
	l1.setPrefix(""); //no prefix
	l1 << " 10." << endl;

	StringList to_validate_list = ListUtils::create<String>(String(stream_by_logger.str()),'\n');
	TEST_EQUAL(to_validate_list.size(),10)

	StringList regex_list;
	regex_list.push_back("DEVELOPMENT  2\\.");
	regex_list.push_back("[0-2][0-9]:[0-5][0-9]:[0-5][0-9]  3\\.");
	regex_list.push_back("[0-2][0-9]:[0-5][0-9]  4\\.");
  regex_list.push_back("[0-9]+/[0-1][0-9]/[0-3][0-9]  5\\.");
	regex_list.push_back("[0-1][0-9]/[0-3][0-9]  6\\.");
  regex_list.push_back("[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]  7\\.");
	regex_list.push_back("[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]  8\\.");
	regex_list.push_back("%  9\\.");
	regex_list.push_back(" 10\\.");

	for (Size i=0;i<regex_list.size();++i)
  {
    boost::regex rx(regex_list[i].c_str());
    TEST_EQUAL(regex_match(to_validate_list[i], rx), true)
	}
}
END_SECTION

START_SECTION((void setPrefix(const std::ostream &s, const std::string &prefix)))
{
  LogStream l1(new LogStreamBuf());
  ostringstream stream_by_logger;
	ostringstream stream_by_logger_otherprefix;
  l1.insert(stream_by_logger);
  l1.insert(stream_by_logger_otherprefix);
  l1.setPrefix(stream_by_logger_otherprefix, "BLABLA"); //message type ("Error", "Warning", "Information", "-")
  l1.setLevel("DEVELOPMENT");
  l1.setPrefix(stream_by_logger, "%y"); //message type ("Error", "Warning", "Information", "-")
  l1 << "  2." << endl;
  l1.setPrefix(stream_by_logger, "%T"); //time (HH:MM:SS)
  l1 << "  3." << endl;
  l1.setPrefix(stream_by_logger, "%t"); //time in short format (HH:MM)
  l1 << "  4." << endl;
  l1.setPrefix(stream_by_logger, "%D"); //date (YYYY/MM/DD)
  l1 << "  5." << endl;
  l1.setPrefix(stream_by_logger, "%d"); // date in short format (MM/DD)
  l1 << "  6." << endl;
  l1.setPrefix(stream_by_logger, "%S"); //time and date (YYYY/MM/DD, HH:MM:SS)
  l1 << "  7." << endl;
  l1.setPrefix(stream_by_logger, "%s"); //time and date in short format (MM/DD, HH:MM)
  l1 << "  8." << endl;
  l1.setPrefix(stream_by_logger, "%%"); //percent sign (escape sequence)
  l1 << "  9." << endl;
  l1.setPrefix(stream_by_logger, ""); //no prefix
  l1 << " 10." << endl;

	StringList to_validate_list = ListUtils::create<String>(String(stream_by_logger.str()),'\n');
	TEST_EQUAL(to_validate_list.size(),10)
	StringList to_validate_list2 = ListUtils::create<String>(String(stream_by_logger_otherprefix.str()),'\n');
	TEST_EQUAL(to_validate_list2.size(),10)

	StringList regex_list;
	regex_list.push_back("DEVELOPMENT  2\\.");
  regex_list.push_back("[0-2][0-9]:[0-5][0-9]:[0-5][0-9]  3\\.");
  regex_list.push_back("[0-2][0-9]:[0-5][0-9]  4\\.");
  regex_list.push_back("[0-9]+/[0-1][0-9]/[0-3][0-9]  5\\.");
  regex_list.push_back("[0-1][0-9]/[0-3][0-9]  6\\.");
  regex_list.push_back("[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]  7\\.");
  regex_list.push_back("[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]  8\\.");
	regex_list.push_back("%  9\\.");
	regex_list.push_back(" 10\\.");

	String other_stream_regex = "BLABLA [ 1][0-9]\\.";
  boost::regex rx2(other_stream_regex);
  // QRegExp rx2(other_stream_regex.c_str());
  // QRegExpValidator v2(rx2, 0);

	for (Size i=0;i<regex_list.size();++i)
	{
    boost::regex rx(regex_list[i].c_str());
    TEST_EQUAL(regex_match(to_validate_list[i], rx), true)
    TEST_EQUAL(regex_match(to_validate_list2[i], rx2), true)
	}

}
END_SECTION

START_SECTION((void flush()))
{
	LogStream l1(new LogStreamBuf());
	ostringstream stream_by_logger;
	l1.insert(stream_by_logger);
	l1 << "flushtest" << endl;
	TEST_EQUAL(stream_by_logger.str(),"flushtest\n")
	l1 << "unfinishedline...\n";
	TEST_EQUAL(stream_by_logger.str(),"flushtest\n")
	l1.flush();
	TEST_EQUAL(stream_by_logger.str(),"flushtest\nunfinishedline...\n")

}
END_SECTION

START_SECTION(([EXTRA]Test minimum string length of output))
{
  // taken from BALL tests, it seems that it checks if the logger crashs if one
  // uses longer lines
  NOT_TESTABLE
  LogStream l1(new LogStreamBuf());
  l1 << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << endl;
}
END_SECTION

START_SECTION(([EXTRA]Test log caching))
{
  String filename;
  NEW_TMP_FILE(filename)
  ofstream s(filename.c_str(), std::ios::out);
  {
    LogStream l1(new LogStreamBuf());
    l1.insert(s);

    l1 << "This is a repeptitive message" << endl;
    l1 << "This is another repeptitive message" << endl;
    l1 << "This is a repeptitive message" << endl;
    l1 << "This is another repeptitive message" << endl;
    l1 << "This is a repeptitive message" << endl;
    l1 << "This is another repeptitive message" << endl;
    l1 << "This is a non-repetitive message" << endl;
  }

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_caching.txt"))
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_FATAL_ERROR))
{
  // remove cout/cerr streams from global instances
  // and append trackable ones
  OpenMS_Log_fatal.remove(cerr);
  ostringstream stream_by_logger;
  {
    OpenMS_Log_fatal.insert(stream_by_logger);

    OPENMS_LOG_FATAL_ERROR << "1\n";
    OPENMS_LOG_FATAL_ERROR << "2" << endl;
  }

  StringList to_validate_list = ListUtils::create<String>(String(stream_by_logger.str()),'\n');
  TEST_EQUAL(to_validate_list.size(),3)

  boost::regex rx(R"(.*LogStream_test\.cpp\(\d+\): \d)");
  for (Size i=0;i<to_validate_list.size() - 1;++i) // there is an extra line since we ended with endl
  {
    TEST_TRUE(regex_search(to_validate_list[i], rx))
  }
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_ERROR))
{
  // remove cout/cerr streams from global instances
  // and append trackable ones
  OpenMS_Log_error.remove(cerr);
  String filename;
  NEW_TMP_FILE(filename)
  ofstream s(filename.c_str(), std::ios::out);
  {
    OpenMS_Log_error.insert(s);

    OPENMS_LOG_ERROR << "1\n";
    OPENMS_LOG_ERROR << "2" << endl;
  }
  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general_red.txt"))
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_WARN))
{
  // remove cout/cerr streams from global instances
  // and append trackable ones
  OpenMS_Log_warn.remove(cout);
  String filename;
  NEW_TMP_FILE(filename)
  ofstream s(filename.c_str(), std::ios::out);
  {
    OpenMS_Log_warn.insert(s);

    OPENMS_LOG_WARN << "1\n";
    OPENMS_LOG_WARN << "2" << endl;
  }
  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general_yellow.txt"))
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_INFO))
{
  // remove cout/cerr streams from global instances
  // and append trackable ones
  OpenMS_Log_info.remove(cout);

  // clear cache to avoid pollution of the test output
  // by previous tests
  OpenMS_Log_info.rdbuf()->clearCache();

  String filename;
  NEW_TMP_FILE(filename)
  ofstream s(filename.c_str(), std::ios::out);
  {
    OpenMS_Log_info.insert(s);

    OPENMS_LOG_INFO << "1\n";
    OPENMS_LOG_INFO << "2" << endl;
  }
  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_DEBUG))
{
  // remove cout/cerr streams from global instances
  // and append trackable ones
  OpenMS_Log_debug.remove(cout);

  // clear cache to avoid pollution of the test output
  // by previous tests
  OpenMS_Log_debug.rdbuf()->clearCache();

  ostringstream stream_by_logger;
  {
    OpenMS_Log_debug.insert(stream_by_logger);

    OPENMS_LOG_DEBUG << "1\n";
    OPENMS_LOG_DEBUG << "2" << endl;
  }

  StringList to_validate_list = ListUtils::create<String>(String(stream_by_logger.str()),'\n');
  TEST_EQUAL(to_validate_list.size(), 3)

  boost::regex rx(R"(.*LogStream_test\.cpp\(\d+\): \d)");
  for (Size i=0;i<to_validate_list.size() - 1;++i) // there is an extra line since we ended with endl
  {
    TEST_TRUE(regex_search(to_validate_list[i], rx))
  }
}
END_SECTION

START_SECTION(([EXTRA] Test caching of empty lines))
{
  ostringstream stream_by_logger;
  {
		LogStream l1(new LogStreamBuf());
		l1.insert(stream_by_logger);
		l1 << "No caching for the following empty lines" << std::endl;
		l1 << "\n\n\n" << std::endl;
	}
	TEST_EQUAL(stream_by_logger.str(), "No caching for the following empty lines\n\n\n\n\n")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



