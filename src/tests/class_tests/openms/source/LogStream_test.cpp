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
#include <iostream>
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
  // Test thread-local logging with OpenMP.
  // Note: Thread-local streams COPY the stream_list_ from global at initialization,
  // so we can't easily add a capture stream and expect thread-local loggers to use it.
  // Instead, we just verify that parallel logging doesn't crash or corrupt data.

  // Test 1: Basic parallel logging to cout (default stream)
  {
    const int num_iterations = 100;

    #ifdef _OPENMP
    omp_set_num_threads(4);
    #pragma omp parallel for
    #endif
    for (int i = 0; i < num_iterations; ++i)
    {
      // Each thread uses its own thread-local LogStream
      OPENMS_LOG_INFO << "iteration_" << i << endl;
    }
    // If we get here without crashing, the test passes
    TEST_EQUAL(true, true)
  }

  // Test 2: High-volume logging stress test
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
    // If we get here without crashing, the test passes
    TEST_EQUAL(true, true)
  }

  // Test 3 (issue #9515): concurrent WARN logging must not crash or tear lines.
  // This faithfully exercises the path that corrupted the heap on Windows (0xC0000374):
  // OPENMS_LOG_WARN writes to the shared std::cerr sink through the shared 'yellow' Colorizer
  // and (pre-fix) through a shared 'static char buf' in LogStreamBuf::syncLF_(), from inside an
  // OpenMP parallel region. We redirect std::cerr's buffer to a temp file (the WARN sink OBJECT
  // is unchanged, so every thread-local buffer still writes to it) to capture output and avoid
  // console spam, then verify each captured line is an intact, well-formed message.
  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename)

    const int N = 20000;
    {
      std::ofstream capture(tmp_filename.c_str());
      std::streambuf* old_cerr = std::cerr.rdbuf(capture.rdbuf());

      #ifdef _OPENMP
      omp_set_num_threads(8);
      #pragma omp parallel for
      #endif
      for (int i = 0; i < N; ++i)
      {
        // Unique-per-iteration: defeats the 2-slot dedup cache so every emit actually reaches
        // distribute_() (the corrupting path). std::endl flushes each message immediately, so the
        // per-thread LogStreamBuf retains no un-flushed tail when we restore std::cerr below
        // (a tail would otherwise be flushed at thread exit, after the restore, and leak to cerr).
        OPENMS_LOG_WARN << "warn_thread_msg_" << i << std::endl;
      }

      std::cerr.flush();
      std::cerr.rdbuf(old_cerr); // restore BEFORE 'capture' is destroyed
    } // 'capture' flushed & closed here

    // Reaching here at all is the primary regression signal (pre-fix: Windows heap abort).
    TEST_EQUAL(true, true)

    // Integrity: every non-empty captured line, with any ANSI color codes stripped, must be a
    // complete "warn_thread_msg_<n>" message. A torn/interleaved line fails the match.
    std::ifstream in(tmp_filename.c_str());
    std::string line;
    Size good = 0;
    boost::regex ansi("\033\\[[0-9;]*m");        // strip color escapes (emitted only on a TTY)
    boost::regex msg("^warn_thread_msg_[0-9]+$");
    while (std::getline(in, line))
    {
      if (line.empty()) { continue; }
      std::string clean = boost::regex_replace(line, ansi, "");
      TEST_EQUAL(boost::regex_match(clean, msg), true)
      ++good;
    }
    // Exactly N lines must survive: messages are unique (no dedup), the default WARN config has a
    // single sink, and std::endl flushes each line, so capture is complete and deterministic.
    // good < N would mean silent message loss (a corruption mode short of a torn line).
    TEST_EQUAL(good == (Size)N, true)
  }
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
  std::string filename;
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

START_SECTION(([EXTRA] LogSinkGuard - RAII removal and re-insertion))
{
  // Test 1: Normal scope exit - guard removes on construction, re-inserts on destruction
  {
    LogStream l1(new LogStreamBuf());
    ostringstream s;
    l1.insert(s);
    l1 << "before_guard" << endl;
    TEST_EQUAL(s.str(), "before_guard\n")

    {
      LogSinkGuard guard(l1, s); // guard removes s immediately
      l1 << "while_guarded" << endl;
      TEST_EQUAL(s.str(), "before_guard\n") // no change, stream was removed by guard
    } // guard destructor re-inserts s

    l1 << "after_guard" << endl;
    TEST_EQUAL(s.str(), "before_guard\nafter_guard\n") // stream is back
  }

  // Test 2: Exception safety - stream should be re-inserted even on exception
  {
    LogStream l1(new LogStreamBuf());
    ostringstream s;
    l1.insert(s);

    try
    {
      LogSinkGuard guard(l1, s); // guard removes s
      l1 << "in_try" << endl;
      throw std::runtime_error("test exception");
    }
    catch (const std::exception&)
    {
      // guard destructor should have run, re-inserting s
    }

    l1 << "after_exception" << endl;
    TEST_EQUAL(s.str(), "after_exception\n") // stream was re-inserted despite exception
  }

  // Test 3: Multiple guards on same stream (nested removal/insertion)
  {
    LogStream l1(new LogStreamBuf());
    ostringstream s;
    l1.insert(s);

    {
      LogSinkGuard guard1(l1, s); // guard1 removes s
      {
        LogSinkGuard guard2(l1, s); // guard2 removes s (already removed - no-op)
        l1 << "deeply_removed" << endl;
      } // guard2 re-inserts s
      l1 << "once_reinserted" << endl;
    } // guard1 re-inserts s (already present - safe/idempotent)
    l1 << "final" << endl;
    TEST_EQUAL(s.str(), "once_reinserted\nfinal\n")
  }
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

	StringList to_validate_list = ListUtils::create<std::string>(stream_by_logger.str(),'\n');
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

	StringList to_validate_list = ListUtils::create<std::string>(stream_by_logger.str(),'\n');
	TEST_EQUAL(to_validate_list.size(),10)
	StringList to_validate_list2 = ListUtils::create<std::string>(stream_by_logger_otherprefix.str(),'\n');
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

	std::string other_stream_regex = "BLABLA [ 1][0-9]\\.";
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
  std::string filename;
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
  // remove cout/cerr streams from the appropriate logger
  // and append trackable ones
  // NOTE: clearCache() outputs cached messages, so call it BEFORE inserting test stream
  ostringstream stream_by_logger;
  {
    getThreadLocalLogFatal().rdbuf()->clearCache();  // outputs to old streams, then clears
    getThreadLocalLogFatal().removeAllStreams();
    getThreadLocalLogFatal().insert(stream_by_logger);

    OPENMS_LOG_FATAL_ERROR << "1\n";
    OPENMS_LOG_FATAL_ERROR << "2" << endl;

    getThreadLocalLogFatal().remove(stream_by_logger);
  }

  StringList to_validate_list = ListUtils::create<std::string>(stream_by_logger.str(),'\n');
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
  // remove cout/cerr streams from the appropriate logger
  // and append trackable ones
  // NOTE: clearCache() outputs cached messages, so call it BEFORE inserting test stream
  std::string filename;
  NEW_TMP_FILE(filename)
  ofstream s(filename.c_str(), std::ios::out);
  {
    getThreadLocalLogError().rdbuf()->clearCache();  // outputs to old streams, then clears
    getThreadLocalLogError().removeAllStreams();
    getThreadLocalLogError().insert(s);

    OPENMS_LOG_ERROR << "1\n";
    OPENMS_LOG_ERROR << "2" << endl;

    getThreadLocalLogError().remove(s);
  }
  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general_red.txt"))
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_WARN))
{
  // remove cout/cerr streams from the appropriate logger
  // and append trackable ones
  // NOTE: clearCache() outputs cached messages, so call it BEFORE inserting test stream
  std::string filename;
  NEW_TMP_FILE(filename)
  ofstream s(filename.c_str(), std::ios::out);
  {
    getThreadLocalLogWarn().rdbuf()->clearCache();  // outputs to old streams, then clears
    getThreadLocalLogWarn().removeAllStreams();
    getThreadLocalLogWarn().insert(s);

    OPENMS_LOG_WARN << "1\n";
    OPENMS_LOG_WARN << "2" << endl;

    getThreadLocalLogWarn().remove(s);
  }
  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general_yellow.txt"))
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_INFO))
{
  // remove cout/cerr streams from the appropriate logger
  // and append trackable ones
  // NOTE: clearCache() outputs cached messages, so call it BEFORE inserting test stream
  std::string filename;
  NEW_TMP_FILE(filename)
  ofstream s(filename.c_str(), std::ios::out);
  {
    getThreadLocalLogInfo().rdbuf()->clearCache();  // outputs to old streams, then clears
    getThreadLocalLogInfo().removeAllStreams();
    getThreadLocalLogInfo().insert(s);

    OPENMS_LOG_INFO << "1\n";
    OPENMS_LOG_INFO << "2" << endl;

    getThreadLocalLogInfo().remove(s);
  }
  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION(([EXTRA] Macro test - OPENMS_LOG_DEBUG))
{
  // remove cout/cerr streams from the appropriate logger
  // and append trackable ones
  // NOTE: clearCache() outputs cached messages, so call it BEFORE inserting test stream
  ostringstream stream_by_logger;
  {
    getThreadLocalLogDebug().rdbuf()->clearCache();  // outputs to old streams, then clears
    getThreadLocalLogDebug().removeAllStreams();
    getThreadLocalLogDebug().insert(stream_by_logger);

    OPENMS_LOG_DEBUG << "1\n";
    OPENMS_LOG_DEBUG << "2" << endl;

    getThreadLocalLogDebug().remove(stream_by_logger);
  }

  StringList to_validate_list = ListUtils::create<std::string>(stream_by_logger.str(),'\n');
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



