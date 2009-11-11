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
#include <QRegExpValidator>

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

START_SECTION((void insertNotification(std::ostream &s, LogStreamNotifier &target, LogLevel min_level=LogStreamBuf::MIN_LEVEL, LogLevel max_level=LogStreamBuf::MAX_LEVEL)))
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

START_SECTION((void setPrefix(const std::string &prefix)))
{
	LogStream l1(new LogStreamBuf());
	ostringstream stream_by_logger;
	l1.insert(stream_by_logger);
	l1.setLevel(DEVELOPMENT);
	l1.setPrefix("%y"); //message type ("Error", "Warning", "Information", "-")
	l1 << "  2." << endl;
	l1.setPrefix("%T"); //time (HH:MM:SS)
	l1 << "  3." << endl;
	l1.setPrefix( "%t"); //time in short format (HH:MM)
	l1 << "  4." << endl;
	l1.setPrefix("%D"); //date (DD.MM.YYYY)
	l1 << "  5." << endl;
	l1.setPrefix("%d"); // date in short format (DD.MM.)
	l1 << "  6." << endl;
	l1.setPrefix("%S"); //time and date (DD.MM.YYYY, HH:MM:SS)
	l1 << "  7." << endl;
	l1.setPrefix("%s"); //time and date in short format (DD.MM., HH:MM)
	l1 << "  8." << endl;
	l1.setPrefix("%%"); //percent sign (escape sequence)
	l1 << "  9." << endl;
	l1.setPrefix(""); //no prefix
	l1 << " 10." << endl;

	StringList to_validate_list = StringList::create(String(stream_by_logger.str()),'\n');
	TEST_EQUAL(to_validate_list.size(),10)

	StringList regex_list;
	regex_list.push_back("DEVELOPMENT  2\\.");
	regex_list.push_back("[0-2][0-9]:[0-5][0-9]:[0-5][0-9]  3\\.");
	regex_list.push_back("[0-2][0-9]:[0-5][0-9]  4\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\.[0-9]+  5\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\.  6\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\.[0-9]+, [0-2][0-9]:[0-5][0-9]:[0-5][0-9]  7\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\., [0-2][0-9]:[0-5][0-9]  8\\.");
	regex_list.push_back("%  9\\.");
	regex_list.push_back(" 10\\.");

	int pos(0);
	for (Size i=0;i<regex_list.size();++i)
		{
		QRegExp rx(regex_list[i].c_str());
		QRegExpValidator v(rx, 0);
		QString to_validate = to_validate_list[i].toQString();
		TEST_EQUAL(v.validate(to_validate,pos)==QValidator::Acceptable, true)
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
  l1.setLevel(DEVELOPMENT);
  l1.setPrefix(stream_by_logger, "%y"); //message type ("Error", "Warning", "Information", "-")
  l1 << "  2." << endl;
  l1.setPrefix(stream_by_logger, "%T"); //time (HH:MM:SS)
  l1 << "  3." << endl;
  l1.setPrefix(stream_by_logger, "%t"); //time in short format (HH:MM)
  l1 << "  4." << endl;
  l1.setPrefix(stream_by_logger, "%D"); //date (DD.MM.YYYY)
  l1 << "  5." << endl;
  l1.setPrefix(stream_by_logger, "%d"); // date in short format (DD.MM.)
  l1 << "  6." << endl;
  l1.setPrefix(stream_by_logger, "%S"); //time and date (DD.MM.YYYY, HH:MM:SS)
  l1 << "  7." << endl;
  l1.setPrefix(stream_by_logger, "%s"); //time and date in short format (DD.MM., HH:MM)
  l1 << "  8." << endl;
  l1.setPrefix(stream_by_logger, "%%"); //percent sign (escape sequence)
  l1 << "  9." << endl;
  l1.setPrefix(stream_by_logger, ""); //no prefix
  l1 << " 10." << endl;
	
	StringList to_validate_list = StringList::create(String(stream_by_logger.str()),'\n');
	TEST_EQUAL(to_validate_list.size(),10)
	StringList to_validate_list2 = StringList::create(String(stream_by_logger_otherprefix.str()),'\n');
	TEST_EQUAL(to_validate_list2.size(),10)

	StringList regex_list;
	regex_list.push_back("DEVELOPMENT  2\\.");
	regex_list.push_back("[0-2][0-9]:[0-5][0-9]:[0-5][0-9]  3\\.");
	regex_list.push_back("[0-2][0-9]:[0-5][0-9]  4\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\.[0-9]+  5\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\.  6\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\.[0-9]+, [0-2][0-9]:[0-5][0-9]:[0-5][0-9]  7\\.");
	regex_list.push_back("[0-3][0-9]\\.[0-1][0-9]\\., [0-2][0-9]:[0-5][0-9]  8\\.");
	regex_list.push_back("%  9\\.");
	regex_list.push_back(" 10\\.");
	
	String other_stream_regex = "BLABLA [ 1][0-9]\\.";
	QRegExp rx2(other_stream_regex.c_str());
	QRegExpValidator v2(rx2, 0);
	
	int pos(0);
	for (Size i=0;i<regex_list.size();++i)
	{
		QRegExp rx(regex_list[i].c_str());
		QRegExpValidator v(rx, 0);
		QString to_validate = to_validate_list[i].toQString();
		QString to_validate2 = to_validate_list2[i].toQString();
		TEST_EQUAL(v.validate(to_validate,pos)==QValidator::Acceptable, true)
		TEST_EQUAL(v2.validate(to_validate2,pos)==QValidator::Acceptable, true)
		
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

START_SECTION(([EXTRA] String LogLevelToStringUpper(LogLevel level)))
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

START_SECTION(([EXTRA] String LogLevelToString(LogLevel level)))
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



