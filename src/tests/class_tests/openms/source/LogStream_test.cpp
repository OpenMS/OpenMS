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
#include <OpenMS/CONCEPT/Colorizer.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <condition_variable>
#include <fstream>
#include <mutex>
#include <sstream>
#include <thread>
#include <boost/regex.hpp>

// OpenMP support
#ifdef _OPENMP
	#include <omp.h>
#endif


///////////////////////////

using namespace OpenMS;
using namespace Logger;
using namespace std;

namespace
{
  /// The text of a destination of a colored stream (debug, warning), without the color codes
  std::string withoutColor(const ostringstream& s)
  {
    return boost::regex_replace(s.str(), boost::regex("\x1b\\[[0-9;]*m"), "");
  }

  /// Logs from a static destructor, after the main thread's thread_local log streams are destroyed (crashed before)
  struct LogsAtExit
  {
    ~LogsAtExit() { OPENMS_LOG_INFO << "LogStream_test: message from a static destructor at exit" << std::endl; }
  } logs_at_exit;
}

START_TEST(LogStream, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] OpenMP - test))
{
  // Test thread-local logging with OpenMP.
  // We just verify that parallel logging doesn't crash or corrupt data.

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
}
END_SECTION

START_SECTION(([EXTRA] thread-safe logging from an OpenMP parallel region (issue #9515)))
{
  // Regression test for https://github.com/OpenMS/OpenMS/issues/9515:
  // emitting warnings from inside an OpenMP parallel region must not corrupt the
  // heap. Before the fix, LogStreamBuf::distribute_()/syncLF_() concurrently
  //  (a) raced on a function-local 'static' line-assembly buffer,
  //  (b) wrote to the shared sink (std::cerr) without synchronization, and
  //  (c) mutated the shared global Colorizer 'yellow' from multiple threads.
  // We capture cerr (the default WARN sink, with the 'yellow' colorizer) by
  // swapping its rdbuf, then hammer it with unique messages from many threads.
  // With the fix the writes are serialized and every message arrives exactly
  // once and intact; without it the run corrupts the heap / interleaves output.

  // make sure this thread's WARN logger writes to cerr (other sections may have
  // reconfigured it) and start from a clean cache
  getThreadLocalLogWarn().rdbuf()->clearCache();
  getThreadLocalLogWarn().insert(std::cerr); // idempotent if already present

  // redirect cerr into a capture buffer (the std::cerr object - and thus every
  // thread-local WARN buffer that points at it - keeps writing there)
  std::ostringstream capture;
  std::streambuf* old_cerr = std::cerr.rdbuf(capture.rdbuf());

  const int num_iterations = 5000;
  #ifdef _OPENMP
  omp_set_num_threads(8);
  #pragma omp parallel for
  #endif
  for (int i = 0; i < num_iterations; ++i)
  {
    OPENMS_LOG_WARN << "racing_line_" << i << std::endl;
  }

  std::cerr.rdbuf(old_cerr); // restore before any assertion/output

  // every unique message must have been distributed exactly once and intact.
  // A plain occurrence count could be fooled by one dropped + one duplicated
  // message cancelling out, so verify each id 0..num_iterations-1 appears once.
  const std::string out = capture.str();
  std::vector<int> seen((Size)num_iterations, 0);
  Size out_of_range = 0;
  boost::regex rx("racing_line_([0-9]+)");
  for (boost::sregex_iterator it(out.begin(), out.end(), rx), rx_end; it != rx_end; ++it)
  {
    const int id = std::stoi((*it)[1].str());
    if (id >= 0 && id < num_iterations) { ++seen[(Size)id]; }
    else { ++out_of_range; }
  }
  Size missing = 0, duplicated = 0;
  for (int v : seen)
  {
    if (v == 0) { ++missing; }
    else if (v > 1) { ++duplicated; }
  }
  TEST_EQUAL(out_of_range, 0)
  TEST_EQUAL(missing, 0)
  TEST_EQUAL(duplicated, 0)
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
        LogSinkGuard guard2(l1, s); // s is already gone, so guard2 has nothing to guard
        l1 << "deeply_removed" << endl;
      } // guard2 must NOT re-insert s -- it never removed it, and guard1 is still active
      l1 << "still_suppressed_by_guard1" << endl;
    } // guard1 re-inserts s
    l1 << "final" << endl;
    // Previously this read "once_reinserted\nfinal\n": the inner guard re-attached the sink and
    // the rest of the outer scope leaked to it. That expectation encoded the bug, not the contract.
    TEST_EQUAL(s.str(), "final\n")
  }

  // Test 4: a message that is never flushed by the writer. This is how OpenMS actually logs --
  // OPENMS_LOG_* messages end in '\n', not std::endl -- so the text is still in the buffer when
  // the guard goes out of scope. It must be discarded there, not handed to the next flush.
  {
    LogStream l1(new LogStreamBuf());
    ostringstream s;
    l1.insert(s);

    {
      LogSinkGuard guard(l1, s);
      l1 << "unflushed_while_guarded\n"; // no endl: nothing is written yet
    } // guard re-inserts s -- but only after dropping the pending text

    l1 << "after_guard" << endl;
    TEST_EQUAL(s.str(), "after_guard\n") // the guarded message must not resurface here
  }

  // Test 5: suppressing one sink must not suppress the others -- a message logged while cout is
  // guarded is not "cancelled", it simply does not go to cout.
  {
    LogStream l1(new LogStreamBuf());
    ostringstream guarded, kept;
    l1.insert(guarded);
    l1.insert(kept);

    {
      LogSinkGuard guard(l1, guarded);
      l1 << "one_sink_only\n";
    }

    TEST_EQUAL(guarded.str(), "")
    TEST_EQUAL(kept.str(), "one_sink_only\n")
  }

  // Test 6: a message with no trailing newline, with a second sink attached. Not every OpenMS log
  // statement terminates its line (e.g. RWrapper's "Running R script ..."), and an unterminated
  // message is parked in the buffer's incomplete_line_ rather than distributed. If the guard only
  // flushed complete lines, that text would outlive it and prefix the next line the restored sink
  // receives.
  {
    LogStream l1(new LogStreamBuf());
    ostringstream guarded, kept;
    l1.insert(guarded);
    l1.insert(kept);

    {
      LogSinkGuard guard(l1, guarded);
      l1 << "unterminated_while_guarded"; // no '\n' at all
    }

    l1 << "after_guard\n";
    l1.flush();
    TEST_EQUAL(guarded.str(), "after_guard\n")            // not "unterminated_while_guardedafter_guard\n"
    TEST_EQUAL(kept.str(), "unterminated_while_guarded\nafter_guard\n") // never guarded, still receives it
  }

  // Test 7: output pending from BEFORE the guard belongs to the guarded sink and must not be
  // swallowed by the guard's clean-up. The buffer is shared by all sinks, so it is delivered when
  // the guard is entered (hence a line of its own) rather than merged with what the scope logs.
  {
    LogStream l1(new LogStreamBuf());
    ostringstream guarded;
    l1.insert(guarded);

    l1 << "pending_before_guard"; // no '\n' yet
    {
      LogSinkGuard guard(l1, guarded);
      l1 << "suppressed\n";
    }
    l1 << "after_guard\n";
    l1.flush();
    TEST_EQUAL(guarded.str(), "pending_before_guard\nafter_guard\n")
  }

  // Test 8: guarding a sink that is not attached must not attach it. "Temporarily remove" has no
  // meaning for a stream that was never a destination, and adding one would redirect output the
  // caller never asked for (and never take it away again).
  {
    LogStream l1(new LogStreamBuf());
    ostringstream attached, never_attached;
    l1.insert(attached);

    {
      LogSinkGuard guard(l1, never_attached);
      l1 << "while_guarded\n";
    }
    l1 << "after_guard\n";
    l1.flush();
    TEST_EQUAL(never_attached.str(), "")
    TEST_EQUAL(attached.str(), "while_guarded\nafter_guard\n") // unrelated sink unaffected
  }
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

START_SECTION((void setConsoleDebugLogging(bool enabled)))
{
  // This thread's debug stream already exists, so a change to the global stream alone would not reach it.
  setConsoleDebugLogging(false);

  ostringstream captured;
  streambuf* cout_buf = cout.rdbuf(captured.rdbuf());
  setConsoleDebugLogging(true);
  OPENMS_LOG_DEBUG << "main thread enabled" << endl;
  std::thread([] { OPENMS_LOG_DEBUG << "new thread enabled" << endl; }).join();
  setConsoleDebugLogging(false);
  OPENMS_LOG_DEBUG << "main thread disabled" << endl;
  std::thread([] { OPENMS_LOG_DEBUG << "new thread disabled" << endl; }).join();
  cout.rdbuf(cout_buf);

  TEST_TRUE(captured.str().find("main thread enabled") != std::string::npos)
  TEST_TRUE(captured.str().find("new thread enabled") != std::string::npos)
  TEST_TRUE(captured.str().find("main thread disabled") == std::string::npos)
  TEST_TRUE(captured.str().find("new thread disabled") == std::string::npos)
  TEST_FALSE(getGlobalLogDebug().hasStream(cout))
}
END_SECTION

START_SECTION(([EXTRA] thread-local streams follow changes of the global stream))
{
  // The global debug stream has no destinations by default, so these first messages go nowhere.
  OPENMS_LOG_DEBUG_NOFILE << "main thread before the change" << endl;

  // a worker that logged before the change, like a thread of an OpenMP pool
  std::mutex m;
  std::condition_variable cv;
  int step = 0;
  auto advance = [&](int to) { { std::lock_guard<std::mutex> lock(m); step = to; } cv.notify_all(); };
  auto await = [&](int at) { std::unique_lock<std::mutex> lock(m); cv.wait(lock, [&] { return step >= at; }); };
  std::thread worker([&] {
    OPENMS_LOG_DEBUG_NOFILE << "worker before the change" << endl;
    advance(1);
    await(2);
    OPENMS_LOG_DEBUG_NOFILE << "worker after insert" << endl;
    advance(3);
    await(5);
    OPENMS_LOG_DEBUG_NOFILE << "worker after remove" << endl;
  });
  await(1);

  ostringstream captured;
  getGlobalLogDebug().insert(captured);
  OPENMS_LOG_DEBUG_NOFILE << "main thread after insert" << endl;
  advance(2);
  await(3);
  // a worker whose first message comes while the destination is attached
  std::thread late_worker([&] {
    OPENMS_LOG_DEBUG_NOFILE << "late worker after insert" << endl;
    advance(4);
    await(5);
    OPENMS_LOG_DEBUG_NOFILE << "late worker after remove" << endl;
  });
  await(4);
  getGlobalLogDebug().remove(captured);
  OPENMS_LOG_DEBUG_NOFILE << "main thread after remove" << endl;
  advance(5);
  worker.join();
  late_worker.join();

  TEST_TRUE(captured.str().find("main thread after insert") != std::string::npos)
  TEST_TRUE(captured.str().find("worker after insert") != std::string::npos)
  TEST_TRUE(captured.str().find("late worker after insert") != std::string::npos)
  TEST_TRUE(captured.str().find("before the change") == std::string::npos)
  // a removed destination may be destroyed right away, so no thread may write to it anymore
  TEST_TRUE(captured.str().find("after remove") == std::string::npos)
}
END_SECTION

START_SECTION(([EXTRA] changes of a thread-local stream apply to its thread only))
{
  ostringstream global_dest, local_dest, later_dest;
  getGlobalLogDebug().insert(global_dest);
  getThreadLocalLogDebug().remove(global_dest); // hide a global destination from this thread
  getThreadLocalLogDebug().insert(local_dest);  // add an own destination
  OPENMS_LOG_DEBUG_NOFILE << "main thread with local changes" << endl;
  std::thread([] { OPENMS_LOG_DEBUG_NOFILE << "other thread" << endl; }).join();

  getGlobalLogDebug().insert(later_dest); // later global changes still apply
  OPENMS_LOG_DEBUG_NOFILE << "main thread after global insert" << endl;
  {
    LogSinkGuard guard(getThreadLocalLogDebug(), later_dest);
    OPENMS_LOG_DEBUG_NOFILE << "main thread guarded" << endl;
    std::thread([] { OPENMS_LOG_DEBUG_NOFILE << "other thread while guarded" << endl; }).join();
  }
  getThreadLocalLogDebug().insert(global_dest); // undoes the local removal
  OPENMS_LOG_DEBUG_NOFILE << "main thread shown again" << endl;

  getThreadLocalLogDebug().remove(local_dest);
  getGlobalLogDebug().remove(global_dest);
  getGlobalLogDebug().remove(later_dest);
  OPENMS_LOG_DEBUG_NOFILE << "main thread after cleanup" << endl;

  auto has = [](const ostringstream& s, const std::string& text) { return s.str().find(text) != std::string::npos; };
  TEST_FALSE(has(global_dest, "main thread with local changes"))
  TEST_TRUE(has(global_dest, "other thread"))
  TEST_FALSE(has(global_dest, "main thread after global insert"))
  TEST_TRUE(has(global_dest, "main thread shown again"))
  TEST_TRUE(has(local_dest, "main thread with local changes"))
  TEST_FALSE(has(local_dest, "other thread"))
  TEST_TRUE(has(local_dest, "main thread shown again"))
  TEST_TRUE(has(later_dest, "main thread after global insert"))
  TEST_FALSE(has(later_dest, "main thread guarded"))
  TEST_TRUE(has(later_dest, "other thread while guarded"))
  TEST_TRUE(has(later_dest, "main thread shown again"))
  TEST_FALSE(has(global_dest, "after cleanup") || has(local_dest, "after cleanup") || has(later_dest, "after cleanup"))
}
END_SECTION

START_SECTION(([EXTRA] setPrefix on a thread-local stream applies to its own destinations only))
{
  ostringstream global_dest, local_dest;
  getGlobalLogDebug().insert(global_dest);
  getThreadLocalLogDebug().insert(local_dest);
  getThreadLocalLogDebug().setPrefix(global_dest, "ignored "); // a global destination's prefix is set on the global stream
  getThreadLocalLogDebug().setPrefix(local_dest, "local ");
  OPENMS_LOG_DEBUG_NOFILE << "first" << endl;
  getThreadLocalLogDebug().setPrefix("all local ");
  getGlobalLogDebug().setPrefix(global_dest, "global ");
  OPENMS_LOG_DEBUG_NOFILE << "second" << endl;
  std::thread([] { OPENMS_LOG_DEBUG_NOFILE << "other thread" << endl; }).join();
  getThreadLocalLogDebug().remove(local_dest);
  getGlobalLogDebug().remove(global_dest);

  TEST_TRUE(local_dest.str().find("local first") != std::string::npos)
  TEST_TRUE(local_dest.str().find("all local second") != std::string::npos)
  TEST_TRUE(global_dest.str().find("ignored") == std::string::npos)
  TEST_TRUE(global_dest.str().find("all local") == std::string::npos)
  TEST_TRUE(global_dest.str().find("global second") != std::string::npos)
  TEST_TRUE(global_dest.str().find("global other thread") != std::string::npos)
}
END_SECTION

START_SECTION(([EXTRA] without a debug destination, OPENMS_LOG_DEBUG skips the whole message))
{
  int evaluated = 0;
  auto argument = [&evaluated] { ++evaluated; return "argument"; };
  // The global debug stream has no destinations by default.
  OPENMS_LOG_DEBUG << argument() << endl;
  OPENMS_LOG_DEBUG_NOFILE << argument() << endl;
  TEST_EQUAL(evaluated, 0)

  ostringstream captured;
  getGlobalLogDebug().insert(captured);
  OPENMS_LOG_DEBUG << argument() << " with file" << endl;
  OPENMS_LOG_DEBUG_NOFILE << argument() << " without file" << endl;
  TEST_EQUAL(evaluated, 2)
  getGlobalLogDebug().remove(captured);
  OPENMS_LOG_DEBUG << argument() << " after remove" << endl;
  TEST_EQUAL(evaluated, 2)

  TEST_TRUE(captured.str().find("argument with file") != std::string::npos)
  TEST_TRUE(captured.str().find("argument without file") != std::string::npos)
  TEST_TRUE(captured.str().find("after remove") == std::string::npos)

  // an else after the macro belongs to the enclosing if
  bool else_taken = false;
  if (evaluated == 0) OPENMS_LOG_DEBUG << "not reached" << endl; else else_taken = true;
  TEST_TRUE(else_taken)
}
END_SECTION

START_SECTION(([EXTRA] a thread-local stream without destinations is failed, so that it does not format messages))
{
  ostringstream dest;
  bool failed_without = false, good_after_global_insert = false, failed_after_global_remove = false;
  bool good_after_local_insert = false, failed_after_local_remove = false;
  std::thread([&] {
    Logger::LogStream& warn = getThreadLocalLogWarn();
    warn.removeAllStreams(); // for this thread only
    failed_without = warn.bad();
    OPENMS_LOG_WARN << "not shown " << 1.5 << endl;

    getGlobalLogWarn().insert(dest);
    OPENMS_LOG_WARN << "after global insert " << 2.5 << endl; // the accessor applies the global change
    good_after_global_insert = warn.good();
    getGlobalLogWarn().remove(dest);
    OPENMS_LOG_WARN << "after global remove" << endl;
    failed_after_global_remove = warn.bad();

    warn.insert(dest);
    good_after_local_insert = warn.good();
    warn << "after local insert " << 3.5 << endl;
    warn.remove(dest);
    failed_after_local_remove = warn.bad();
  }).join();

  TEST_TRUE(failed_without)
  TEST_TRUE(good_after_global_insert)
  TEST_TRUE(failed_after_global_remove)
  TEST_TRUE(good_after_local_insert)
  TEST_TRUE(failed_after_local_remove)
  TEST_TRUE(dest.str().find("not shown") == std::string::npos)
  TEST_TRUE(dest.str().find("after global insert 2.5") != std::string::npos)
  TEST_TRUE(dest.str().find("after global remove") == std::string::npos)
  TEST_TRUE(dest.str().find("after local insert 3.5") != std::string::npos)
  // global streams keep their state
  TEST_TRUE(getGlobalLogDebug().good())
}
END_SECTION

START_SECTION(([EXTRA] a local insert adds a destination that the global stream removed in the meantime))
{
  ostringstream sink;
  getGlobalLogDebug().insert(sink);
  getThreadLocalLogDebug().remove(sink);
  getGlobalLogDebug().remove(sink);
  getThreadLocalLogDebug().insert(sink);
  TEST_TRUE(getThreadLocalLogDebug().hasStream(sink))
  OPENMS_LOG_DEBUG_NOFILE << "after local insert" << endl;
  getThreadLocalLogDebug().remove(sink);
  TEST_TRUE(sink.str().find("after local insert") != std::string::npos)
}
END_SECTION

START_SECTION(([EXTRA] repeated messages after a change of the destinations))
{
  // the debug stream is colored
  auto plain = [](const ostringstream& s) { return boost::regex_replace(s.str(), boost::regex("\x1b\\[[0-9;]*m"), ""); };

  // A destination attached after a message gets its repetitions; it gets no repeat count of messages it did not
  // get. A remaining destination gets its pending repeat count. A removed one does not: it may be destroyed already.
  ostringstream first, second, kept;
  getGlobalLogDebug().insert(first);
  getGlobalLogDebug().insert(kept);
  OPENMS_LOG_DEBUG_NOFILE << "same message" << endl;
  OPENMS_LOG_DEBUG_NOFILE << "same message" << endl;
  getGlobalLogDebug().remove(first);
  getGlobalLogDebug().insert(second);
  OPENMS_LOG_DEBUG_NOFILE << "same message" << endl;
  getGlobalLogDebug().remove(second);
  getGlobalLogDebug().remove(kept);
  TEST_EQUAL(plain(first), "same message\n")
  TEST_EQUAL(plain(second), "same message\n")
  TEST_EQUAL(plain(kept), "same message\n<same message> occurred 2 times\nsame message\n")

  // changes of the thread-local stream: pending repeat counts go to the destinations before the change
  ostringstream local, later;
  getThreadLocalLogDebug().insert(local);
  OPENMS_LOG_DEBUG_NOFILE << "local message" << endl;
  OPENMS_LOG_DEBUG_NOFILE << "local message" << endl;
  getThreadLocalLogDebug().insert(later);
  OPENMS_LOG_DEBUG_NOFILE << "local message" << endl;
  getThreadLocalLogDebug().remove(later);
  getThreadLocalLogDebug().remove(local);
  TEST_EQUAL(plain(local), "local message\n<local message> occurred 2 times\nlocal message\n")
  TEST_EQUAL(plain(later), "local message\n")
}
END_SECTION

START_SECTION(([EXTRA] LogSinkGuard restores only what it removed))
{
  // A global destination removed while the guard hid it on this thread stays removed (it may be destroyed).
  ostringstream console;
  getGlobalLogDebug().insert(console);
  {
    LogSinkGuard quiet(getThreadLocalLogDebug(), console);
    getGlobalLogDebug().remove(console);
  }
  OPENMS_LOG_DEBUG_NOFILE << "after the guard" << endl;
  TEST_FALSE(getThreadLocalLogDebug().hasStream(console))
  TEST_TRUE(console.str().find("after the guard") == std::string::npos)

  // A destination of the thread-local stream itself is inserted again, with its prefix.
  ostringstream own;
  getThreadLocalLogDebug().insert(own);
  getThreadLocalLogDebug().setPrefix(own, "own: ");
  {
    LogSinkGuard quiet(getThreadLocalLogDebug(), own);
    OPENMS_LOG_DEBUG_NOFILE << "guarded" << endl;
  }
  OPENMS_LOG_DEBUG_NOFILE << "restored" << endl;
  getThreadLocalLogDebug().remove(own);
  TEST_EQUAL(withoutColor(own), "own: restored\n")
}
END_SECTION

START_SECTION(([EXTRA] a LogSinkGuard on a thread-local stream suppresses the stream for its whole scope))
{
  ostringstream sink;
  getGlobalLogDebug().insert(sink);
  {
    LogSinkGuard outer(getThreadLocalLogDebug(), sink);
    getGlobalLogDebug().remove(sink);
    getGlobalLogDebug().insert(sink); // inserted again: still suppressed on this thread
    OPENMS_LOG_DEBUG_NOFILE << "1: inside, after the global insertion" << endl;
    {
      LogSinkGuard inner(getThreadLocalLogDebug(), sink); // does nothing: suppressed already
    }
    OPENMS_LOG_DEBUG_NOFILE << "2: inside, after an inner guard" << endl;
  }
  OPENMS_LOG_DEBUG_NOFILE << "3: after the guard" << endl;

  // a removal in the guarded scope lasts beyond it
  {
    LogSinkGuard quiet(getThreadLocalLogDebug(), sink);
    getThreadLocalLogDebug().remove(sink);
  }
  OPENMS_LOG_DEBUG_NOFILE << "4: removed in the guarded scope" << endl;
  getThreadLocalLogDebug().insert(sink);
  getGlobalLogDebug().remove(sink);
  TEST_EQUAL(withoutColor(sink), "3: after the guard\n")

  // an own destination inserted again in the guarded scope stays a single destination
  ostringstream own;
  getThreadLocalLogDebug().insert(own);
  {
    LogSinkGuard quiet(getThreadLocalLogDebug(), own);
    getThreadLocalLogDebug().insert(own);
  }
  OPENMS_LOG_DEBUG_NOFILE << "once" << endl;
  getThreadLocalLogDebug().remove(own);
  TEST_EQUAL(withoutColor(own), "once\n")
}
END_SECTION

START_SECTION(([EXTRA] a LogSinkGuard on a global stream keeps a destination hidden on a thread))
{
  ostringstream sink;
  getGlobalLogDebug().insert(sink);
  getThreadLocalLogDebug().remove(sink); // hidden on this thread
  {
    LogSinkGuard quiet(getGlobalLogDebug(), sink); // e.g. in a library function
  }
  OPENMS_LOG_DEBUG_NOFILE << "1: after a global guard" << endl;
  getThreadLocalLogDebug().insert(sink); // shown again
  {
    LogSinkGuard outer(getThreadLocalLogDebug(), sink);
    {
      LogSinkGuard inner(getGlobalLogDebug(), sink);
    }
    OPENMS_LOG_DEBUG_NOFILE << "2: inside the outer guard" << endl;
  }
  OPENMS_LOG_DEBUG_NOFILE << "3: after the outer guard" << endl;
  getGlobalLogDebug().remove(sink);
  TEST_EQUAL(withoutColor(sink), "3: after the outer guard\n")
}
END_SECTION

START_SECTION(([EXTRA] a stream inserted again is a new destination))
{
  // also a new stream at the address of a destroyed one, which this simulates deterministically
  ostringstream s;
  getGlobalLogDebug().insert(s);
  OPENMS_LOG_DEBUG_NOFILE << "x" << endl;
  getGlobalLogDebug().remove(s);
  s.str("");
  getGlobalLogDebug().insert(s);
  OPENMS_LOG_DEBUG_NOFILE << "x" << endl; // not a repeat for the new destination

  getThreadLocalLogDebug().removeAllStreams(); // hides the current destinations from this thread
  getGlobalLogDebug().remove(s);
  getGlobalLogDebug().insert(s);
  OPENMS_LOG_DEBUG_NOFILE << "y" << endl; // the new destination is not hidden
  getGlobalLogDebug().remove(s);
  TEST_EQUAL(withoutColor(s), "x\ny\n")
}
END_SECTION

START_SECTION(([EXTRA] LogSinkGuard on a global stream drains the stream of the calling thread))
{
  ostringstream other, sink;
  getGlobalLogDebug().insert(other);
  getGlobalLogDebug().insert(sink);
  OPENMS_LOG_DEBUG_NOFILE << "before the guard\n"; // "\n" without endl: still buffered
  {
    LogSinkGuard quiet(getGlobalLogDebug(), sink);
    OPENMS_LOG_DEBUG_NOFILE << "inside 1" << endl;
    OPENMS_LOG_DEBUG_NOFILE << "inside 2\n";
  }
  OPENMS_LOG_DEBUG_NOFILE << "after" << endl;
  OPENMS_LOG_DEBUG_NOFILE << "partial before the guard: " << flush;
  {
    LogSinkGuard quiet(getGlobalLogDebug(), sink);
    OPENMS_LOG_DEBUG_NOFILE << "inside 3" << endl;
  }
  OPENMS_LOG_DEBUG_NOFILE << "rep" << endl;
  OPENMS_LOG_DEBUG_NOFILE << "rep" << endl;
  {
    LogSinkGuard quiet(getGlobalLogDebug(), sink);
    OPENMS_LOG_DEBUG_NOFILE << "inside 4" << endl;
  }
  OPENMS_LOG_DEBUG_NOFILE << "after 3" << endl;
  getGlobalLogDebug().remove(sink);
  getGlobalLogDebug().remove(other);
  TEST_EQUAL(withoutColor(sink), "before the guard\nafter\npartial before the guard: \nrep\n<rep> occurred 2 times\nafter 3\n")
  TEST_TRUE(withoutColor(other).find("inside 2") != std::string::npos)
}
END_SECTION

START_SECTION(([EXTRA] a failing destination with exceptions enabled does not stop logging))
{
  struct FailingBuf : std::streambuf
  {
    int overflow(int) override { return traits_type::eof(); }
    std::streamsize xsputn(const char*, std::streamsize) override { return 0; }
  } failing_buf;
  std::ostream failing(&failing_buf);
  failing.exceptions(std::ios::badbit);
  ostringstream good, third;
  getGlobalLogDebug().insert(failing);
  getGlobalLogDebug().insert(good);
  getGlobalLogDebug().insert(third);
  bool threw = false;
  try
  {
    OPENMS_LOG_DEBUG_NOFILE << "a" << endl;
    OPENMS_LOG_DEBUG_NOFILE << "a" << endl;
    getGlobalLogDebug().remove(third);
    OPENMS_LOG_DEBUG_NOFILE << "b" << endl; // the accessor writes the pending count of "a" first
    OPENMS_LOG_DEBUG_NOFILE << "b" << endl;
    {
      LogSinkGuard quiet(getThreadLocalLogDebug(), good); // writes the pending count of "b"
    }
    OPENMS_LOG_DEBUG_NOFILE << "c" << endl;
  }
  catch (...)
  {
    threw = true;
  }
  getGlobalLogDebug().remove(failing);
  getGlobalLogDebug().remove(good);
  TEST_FALSE(threw)
  TEST_EQUAL(withoutColor(good), "a\n<a> occurred 2 times\nb\n<b> occurred 2 times\nc\n")
}
END_SECTION

START_SECTION(([EXTRA] pending lines and repeat counts stay with the destinations before a change))
{
  // insert()
  ostringstream a, b;
  getGlobalLogDebug().insert(a);
  OPENMS_LOG_DEBUG_NOFILE << "x" << endl;
  OPENMS_LOG_DEBUG_NOFILE << "x\n" << "before insert\n"; // not flushed yet
  getThreadLocalLogDebug().insert(b);
  OPENMS_LOG_DEBUG_NOFILE << "after insert" << endl;
  getThreadLocalLogDebug().remove(b);
  getGlobalLogDebug().remove(a);
  TEST_EQUAL(withoutColor(a), "x\nbefore insert\n<x> occurred 2 times\nafter insert\n")
  TEST_EQUAL(withoutColor(b), "after insert\n")

  // setPrefix() through a kept reference applies a pending global change first, so the pending repeat count
  // keeps the prefix of the lines it counts
  ostringstream own, other;
  Logger::LogStream& log = getThreadLocalLogDebug();
  log.insert(own);
  OPENMS_LOG_DEBUG_NOFILE << "p" << endl;
  OPENMS_LOG_DEBUG_NOFILE << "p" << endl;
  getGlobalLogDebug().insert(other);
  log.setPrefix(std::string("NEW: "));
  log.remove(own);
  getGlobalLogDebug().remove(other);
  TEST_EQUAL(withoutColor(own), "p\n<p> occurred 2 times\n")
}
END_SECTION

START_SECTION(([EXTRA] messages from the destructor of a thread_local object reach their destinations))
{
  ostringstream dest;
  getGlobalLogDebug().insert(dest);
  std::thread([] {
    struct LogsInDestructor
    {
      ~LogsInDestructor() { OPENMS_LOG_DEBUG_NOFILE << "from a thread_local destructor" << endl; }
    };
    thread_local LogsInDestructor object; // constructed before the thread's log stream, so destroyed after it
    (void)object;
    OPENMS_LOG_DEBUG_NOFILE << "worker" << endl;
  }).join();
  getGlobalLogDebug().remove(dest);
  TEST_TRUE(dest.str().find("worker") != std::string::npos)
  TEST_TRUE(dest.str().find("from a thread_local destructor") != std::string::npos)
}
END_SECTION

START_SECTION(([EXTRA] the log levels use their own Colorizers, not the public ones))
{
  // The public Colorizers are static objects of another file: static initializers and destructors may log before they
  // are constructed or after they are destroyed. Logging leaves them alone.
  yellow("prepared");
  ostringstream dest;
  getGlobalLogWarn().insert(dest);
  OPENMS_LOG_WARN << "warning" << endl;
  getGlobalLogWarn().remove(dest);
  TEST_EQUAL(withoutColor(dest), "warning\n")
  ostringstream colored;
  colored << yellow;
  TEST_TRUE(colored.str().find("prepared") != std::string::npos)
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



