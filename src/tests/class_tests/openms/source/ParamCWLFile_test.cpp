// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Chris Bielow $
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#include <fstream>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/ParamCWLFile.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

///////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;

START_TEST(ParamCWLFile, "$Id")

START_SECTION((ParamCWLFile::load))
{
  String filename;
  NEW_TMP_FILE(filename)
  Param param;
  param.setValue("test:1:value", 1, "description");

  // Check that FileNotFound is being thrown
  TEST_EXCEPTION(Exception::FileNotFound, ParamCWLFile::load("/does/not/exist/FileDoesNotExist.json", param))

  // Check parsing error is thrown
  std::ofstream ofs(filename.c_str(), std::ios::out);
  ofs << "not a json";
  ofs.close();

  TEST_EXCEPTION(Exception::ParseError, ParamCWLFile::load(filename.c_str(), param))

  // Check all types can be parsed
  /// set all expected params
  param.setValue("test:1:bool1", "false");
  param.setValidStrings("test:1:bool1", {"true", "false"});
  param.setValue("test:1:bool2", "false");
  param.setValidStrings("test:1:bool2", {"true", "false"});
  param.setValue("test:1:bool3", "true");
  param.setValidStrings("test:1:bool3", {"false", "true"});
  param.setValue("test:1:bool4", "true");
  param.setValidStrings("test:1:bool4", {"false", "true"});
  param.setValue("test:1:int", 0);
  param.setValue("test:1:double", 0.);
  param.setValue("test:1:string", "");
  param.setValue("test:1:int_list", std::vector<int>{});
  param.setValue("test:1:double_list", std::vector<double>{});
  param.setValue("test:1:string_list", std::vector<std::string>{});


  // create matching json file
  ofs.open(filename.c_str(), std::ios::out);
  ofs << "{\n"
         "  \"bool1\": true,\n"
         "  \"bool2\": false,\n"
         "  \"bool3\": true,\n"
         "  \"bool4\": false,\n"
         "  \"int\": 5,\n"
         "  \"double\": 6.1,\n"
         "  \"string\": \"Hello OpenMS\",\n"
         "  \"int_list\": [10, 11, 12],\n"
         "  \"double_list\": [13.25, 15.125],\n"
         "  \"string_list\": [\"SeqAn\", \"rocks\"]\n"
         "}\n";
  ofs.close();
  ParamCWLFile::load(filename.c_str(), param);

  TEST_EQUAL(param.getValue("test:1:bool1").toBool(), true);
  TEST_EQUAL(param.getValue("test:1:bool2").toBool(), false);
  TEST_EQUAL(param.getValue("test:1:bool3").toBool(), true);
  TEST_EQUAL(param.getValue("test:1:bool4").toBool(), false);
  TEST_EQUAL(int(param.getValue("test:1:int")), 5);
  TEST_EQUAL(double(param.getValue("test:1:double")), 6.1);
  TEST_STRING_EQUAL(std::string(param.getValue("test:1:string")), "Hello OpenMS");

  std::vector<int> int_list = param.getValue("test:1:int_list").toIntVector();
  TEST_EQUAL(int_list.size(), 3);
  TEST_EQUAL(int_list[0], 10);
  TEST_EQUAL(int_list[1], 11);
  TEST_EQUAL(int_list[2], 12);

  std::vector<double> double_list = param.getValue("test:1:double_list").toDoubleVector();
  TEST_EQUAL(double_list.size(), 2);
  TEST_EQUAL(double_list[0], 13.25);
  TEST_EQUAL(double_list[1], 15.125);

  std::vector<std::string> string_list = param.getValue("test:1:string_list").toStringVector();
  TEST_EQUAL(string_list.size(), 2);
  TEST_EQUAL(string_list[0], "SeqAn");
  TEST_EQUAL(string_list[1], "rocks");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop
