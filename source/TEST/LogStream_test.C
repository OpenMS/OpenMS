// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Stephan Aiche, Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------


/**

  Generously provided by the BALL people, taken from version 1.2

*/

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/LogStream.h>
///////////////////////////

using namespace OpenMS;
using namespace Logger;
using namespace std;

class TestTarget
  :  public LogStreamNotifier
{
  public:
  virtual void logNotify()
  {
    notified = true;
    return;
  }
  bool notified;
};


START_TEST(LogStream, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(LogStream(LogStreamBuf *buf=0, bool delete_buf=true, bool associate_stdio=false))
{
  LogStream* l1 = 0;
  l1 = new LogStream((LogStreamBuf*)0);
  TEST_NOT_EQUAL(l1, 0)
  delete l1;

  LogStream* l2 = 0;
  LogStreamBuf* lb2 = new LogStreamBuf();
  l2 = new LogStream(lb2);
  TEST_NOT_EQUAL(l2, 0)
  delete l2;
}
END_SECTION

START_SECTION((virtual ~LogStream()))
{
  LogStream* l1 = new LogStream((LogStreamBuf*)0);
  delete l1;
}
END_SECTION
// if we add these two tests to the end of the TEST 
// we have correct out put on the command line
START_SECTION((LogStreamBuf* operator->()))
{
  LogStream l1(new LogStreamBuf());
  l1->sync();
}
END_SECTION

START_SECTION((LogStreamBuf* rdbuf()))
{
  LogStream l1(new LogStreamBuf());
  // small workaround since TEST_NOT_EQUAL(l1.rdbuf, 0) would expand to
  // cout << ls.rdbuf()
  // which kills the cout buffer
  TEST_NOT_EQUAL((l1.rdbuf()==0), true)
}
END_SECTION
  
START_SECTION((void setLevel(LogLevel level)))
{
  String filename;
  NEW_TMP_FILE(filename)
  LogStream l1(new LogStreamBuf());
  ofstream s(filename.c_str(), std::ios::out);
  l1.insert(s, DEVELOPMENT , ERROR);

  l1 << "1" << endl;
  l1.setLevel(INFORMATION);
  l1 << "2" << endl;
  l1.setLevel(FATAL_ERROR);
  l1 << "X" << endl;

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION((LogLevel getLevel()))
{
  LogStream l1(new LogStreamBuf());
  TEST_EQUAL(l1.getLevel(), DEVELOPMENT)
  l1.setLevel(FATAL_ERROR);
  TEST_EQUAL(l1.getLevel(), FATAL_ERROR)
}
END_SECTION

START_SECTION((LogStream& level(LogLevel level)))
{
  String filename;
  NEW_TMP_FILE(filename)
  LogStream l1(new LogStreamBuf());
  ofstream s(filename.c_str(), std::ios::out);
  l1.insert(s, DEVELOPMENT , ERROR);

  l1.level(DEVELOPMENT) << "1" <<endl;
  l1.level(ERROR) << "2" <<endl;
  l1.level(FATAL_ERROR) << "X" <<endl;

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION((void insert(std::ostream &s, LogLevel min_level=LogStreamBuf::MIN_LEVEL, LogLevel max_level=LogStreamBuf::MAX_LEVEL)))
{
  String filename;
  NEW_TMP_FILE(filename)
  LogStream l1(new LogStreamBuf());
  ofstream s(filename.c_str(), std::ios::out);
  l1.insert(s, ERROR, ERROR);

  l1.level(WARNING) << "X" << endl;
  l1.level(ERROR) << "1" << endl;
  l1.level(ERROR)  << "2" << endl;
  l1.level(FATAL_ERROR)<< "X" << endl;

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION((void remove(std::ostream &s)))
{
  LogStream l1(new LogStreamBuf());
  ofstream s;
  l1.insert(s);
  l1.remove(s);
  // make sure we can remove it twice without harm
  l1.remove(s);
}
END_SECTION

START_SECTION((void insertNotification(std::ostream &s, LogStreamNotifier &target, int min_level=LogStreamBuf::MIN_LEVEL, int max_level=LogStreamBuf::MAX_LEVEL)))
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

START_SECTION((void setMinLevel(const std::ostream &s, LogLevel min_level)))
{
  String filename;
  NEW_TMP_FILE(filename)
  LogStream l1(new LogStreamBuf());
  ofstream s(filename.c_str(), std::ios::out);
  l1.insert(s, DEVELOPMENT);
  l1.setMinLevel(s, WARNING);
  l1.level(INFORMATION) << "X" << endl;
  l1.level(WARNING) << "1" << endl;
  l1.level(ERROR) << "2" << endl;

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION((void setMaxLevel(const std::ostream &s, LogLevel max_level)))
{
  String filename;
  NEW_TMP_FILE(filename)
  LogStream l1(new LogStreamBuf());
  ofstream s(filename.c_str(), std::ios::out);
  l1.insert(s, DEVELOPMENT);
  l1.setMaxLevel(s, ERROR);
  l1.level(WARNING) << "1" << endl;
  l1.level(ERROR) << "2" << endl;
  l1.level(FATAL_ERROR) << "X" << endl;

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_general.txt"))
}
END_SECTION

START_SECTION((void setPrefix(const std::ostream &s, const String &prefix)))
{
	/*
  String filename;
  NEW_TMP_FILE(filename)
  LogStream l1(new LogStreamBuf());
  ofstream s(filename.c_str(), std::ios::out);;
  l1.insert(s, DEVELOPMENT);
  l1.setPrefix(s, "%l"); //loglevel
  l1.level(1) << "  1." << endl;
  l1.setPrefix(s, "%y"); //message type ("Error", "Warning", "Information", "-")
  l1.level(2) << "  2." << endl;
  l1.setPrefix(s, "%T"); //time (HH:MM:SS)
  l1.level(3) << "  3." << endl;
  l1.setPrefix(s, "%t"); //time in short format (HH:MM)
  l1.level(4) << "  4." << endl;
  l1.setPrefix(s, "%D"); //date (DD.MM.YYYY)
  l1.level(5) << "  5." << endl;
  l1.setPrefix(s, "%d"); // date in short format (DD.MM.)
  l1.level(6) << "  6." << endl;
  l1.setPrefix(s, "%S"); //time and date (DD.MM.YYYY, HH:MM:SS)
  l1.level(7) << "  7." << endl;
  l1.setPrefix(s, "%s"); //time and date in short format (DD.MM., HH:MM)
  l1.level(8) << "  8." << endl;
  l1.setPrefix(s, "%%"); //percent sign (escape sequence)
  l1.level(9) << "  9." << endl;
  l1.setPrefix(s, ""); //no prefix
  l1.level(10) << " 10." << endl;
	*/
  /*
  TEST_EQUAL(l1.getNumberOfLines(), 10)
  */
  //TEST_FILE_REGEXP(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("LogStream_test_setPrefix.txt"))
}
END_SECTION

START_SECTION((void disableOutput()))
{
  // TODO
}
END_SECTION

START_SECTION((void enableOutput()))
{
  // TODO
}
END_SECTION

START_SECTION((bool outputEnabled() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void flush()))
{
  // TODO
}
END_SECTION

START_SECTION(([EXTRA]Test minimum string length of output))
{
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
    l1.insert(s, DEVELOPMENT);

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

START_SECTION((String LogLevelToStringUpper(LogLevel level)))
{
	TEST_STRING_EQUAL(LogLevelToStringUpper(FATAL_ERROR), "FATAL_ERROR")
	TEST_STRING_EQUAL(LogLevelToStringUpper(ERROR), "ERROR")
	TEST_STRING_EQUAL(LogLevelToStringUpper(WARNING), "WARNING")
	TEST_STRING_EQUAL(LogLevelToStringUpper(INFORMATION), "INFORMATION")
	TEST_STRING_EQUAL(LogLevelToStringUpper(DEBUG), "DEBUG")
	TEST_STRING_EQUAL(LogLevelToStringUpper(DEBUG_INTENSE), "DEBUG_INTENSE")
	TEST_STRING_EQUAL(LogLevelToStringUpper(DEVELOPMENT), "DEVELOPMENT")
}
END_SECTION

START_SECTION((String LogLevelToString(LogLevel level)))
{
	TEST_STRING_EQUAL(LogLevelToString(FATAL_ERROR), "fatal_error")
	TEST_STRING_EQUAL(LogLevelToString(ERROR), "error")
	TEST_STRING_EQUAL(LogLevelToString(WARNING), "warning")
	TEST_STRING_EQUAL(LogLevelToString(INFORMATION), "information")
	TEST_STRING_EQUAL(LogLevelToString(DEBUG), "debug")
	TEST_STRING_EQUAL(LogLevelToString(DEBUG_INTENSE), "debug_intense")
	TEST_STRING_EQUAL(LogLevelToString(DEVELOPMENT), "development")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



