// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/LogConfigHandler.h>
#include <QRegExpValidator>
///////////////////////////

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

  Param p = LogConfigHandler::getInstance().parse(settings);
  
  // p should contain a list of the above set commands  
  StringList parsedConfigs = static_cast<StringList>( p.getValue(LogConfigHandler::PARAM_NAME) );

  TEST_EQUAL(parsedConfigs[0] , "DEBUG add cout FILE")
  TEST_EQUAL(parsedConfigs[1] , "DEBUG add a.out FILE")
  TEST_EQUAL(parsedConfigs[2] , "INFO add a.out FILE")
  TEST_EQUAL(parsedConfigs[3] , "FATAL_ERROR add cerr FILE")
  
  StringList settings2;
  settings2.push_back("DEBUG");
  
  TEST_EXCEPTION(Exception::ParseError, LogConfigHandler::getInstance().parse(settings2));  
}
END_SECTION

START_SECTION((void configure(const Param &param)))
{
  StringList settings;
  settings.push_back("INFO add testing_info_warn_stream STRING");
  settings.push_back("WARNING add testing_info_warn_stream STRING");
  settings.push_back("ERROR add only_error_string_stream STRING");
  settings.push_back("INFO remove cout FILE");
  settings.push_back("WARNING remove cout");
  settings.push_back("ERROR remove cerr FILE");
  
  Param p;
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");
  
  LogConfigHandler::getInstance().configure(p);

  LOG_INFO << "1" << endl;
  LOG_INFO << "2" << endl;
  LOG_WARN << "3" << endl;
  LOG_ERROR << "4" << endl;

  settings.clear();
  settings.push_back("WARNING clear");
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");

  LogConfigHandler::getInstance().configure(p);

  // this should go into nowhere
  LOG_WARN << "5" << endl;

  ostringstream& info_warn_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance().getStream("testing_info_warn_stream"));
  String info_warn_stream_content(info_warn_stream.str());
  StringList info_warn_result;
  info_warn_stream_content.trim().split('\n', info_warn_result, true );

  TEST_EQUAL(info_warn_result.size() , 3)

  // check output with regex

  int pos(0);

  int i = 1;
  for(StringList::ConstIterator it = info_warn_result.begin() ; it != info_warn_result.end(); ++it)
  {
    QString pattern("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] ");
    pattern.append(QString::number(i));
    QRegExp rx(pattern);
    QRegExpValidator v(rx, 0);
    QString to_validate = it->toQString();
    TEST_EQUAL(v.validate(to_validate,pos)==QValidator::Acceptable, true)
    ++i;
  }
  ostringstream& error_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance().getStream("only_error_string_stream"));
  String error_stream_content(error_stream.str());
  StringList error_result;
  error_stream_content.trim().split('\n', error_result, true );


  TEST_EQUAL(error_result.size(), 1)

  QString pattern("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] 4");
  QRegExp rx(pattern);
  QRegExpValidator v(rx, 0);
  QString to_validate = error_result[0].toQString();
  TEST_EQUAL(v.validate(to_validate,pos)==QValidator::Acceptable, true)
}
END_SECTION

START_SECTION((ostream& getStream(const String &stream_name)))
{
  StringList settings;
  settings.push_back("INFO add testing_getStream STRING");

  Param p;
  p.setValue(LogConfigHandler::PARAM_NAME, settings, "List of all settings that should be applied to the current Logging Configuration");

  LogConfigHandler::getInstance().configure(p);

  LOG_INFO << "getStream 1" << endl;

  ostringstream& info_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance().getStream("testing_getStream"));
  String info_content(info_stream.str());

  StringList info_result;
  info_content.trim().split('\n', info_result, true );

  TEST_EQUAL(info_result.size() , 1)

  // check if everything landed in the stream we wanted
  int pos(0);

  QString pattern("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] getStream 1");
  QRegExp rx(pattern);
  QRegExpValidator v(rx, 0);
  QString to_validate = info_result[0].toQString();
  TEST_EQUAL(v.validate(to_validate,pos)==QValidator::Acceptable, true)
}
END_SECTION

START_SECTION((static LogConfigHandler& getInstance()))
{
  TEST_NOT_EQUAL(&LogConfigHandler::getInstance(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
