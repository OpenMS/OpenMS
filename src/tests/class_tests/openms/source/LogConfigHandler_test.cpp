// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

  Param p = LogConfigHandler::getInstance().parse(settings);

  // p should contain a list of the above set commands
  StringList parsedConfigs = p.getValue(LogConfigHandler::PARAM_NAME);

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
  String pattern("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] ");
  boost::regex rx(pattern);

  int i = 1;
  for(StringList::const_iterator it = info_warn_result.begin() ; it != info_warn_result.end(); ++it)
  {
    boost::regex rx(pattern + i);
    TEST_EQUAL(regex_match(*it, rx), true)
    ++i;
  }
  ostringstream& error_stream = static_cast<ostringstream&>(LogConfigHandler::getInstance().getStream("only_error_string_stream"));
  String error_stream_content(error_stream.str());
  StringList error_result;
  error_stream_content.trim().split('\n', error_result, true );


  TEST_EQUAL(error_result.size(), 1)

  String pattern2("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] 4");
  boost::regex rx2(pattern2);
  TEST_EQUAL(regex_match(error_result[0], rx2), true)
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
  String pattern("\\[[0-9]+/[0-1][0-9]/[0-3][0-9], [0-2][0-9]:[0-5][0-9]:[0-5][0-9]\\] getStream 1");
  boost::regex rx(pattern);
  TEST_EQUAL(regex_match(info_result[0], rx), true)
}
END_SECTION

LogConfigHandler* nullPointer = nullptr;
START_SECTION((static LogConfigHandler& getInstance()))
{
  TEST_NOT_EQUAL(&LogConfigHandler::getInstance(), nullPointer)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
