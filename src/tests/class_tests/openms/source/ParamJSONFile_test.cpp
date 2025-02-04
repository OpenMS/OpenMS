// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <fstream>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParamJSONFile.h>

///////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;

START_TEST(ParamJSONFile, "$Id")

START_SECTION((bool ParamJSONFile::load(const std::string& filename, Param& param)))
{
  String filename;
  NEW_TMP_FILE(filename)
  Param param;
  param.setValue("test:1:value", 1, "description");

  // Check that FileNotFound is being thrown
  TEST_EXCEPTION(Exception::FileNotFound, ParamJSONFile::load("/does/not/exist/FileDoesNotExist.json", param))

  // Check parsing error is thrown
  std::ofstream ofs(filename.c_str(), std::ios::out);
  ofs << "not a json";
  ofs.close();

  TEST_EXCEPTION(Exception::ParseError, ParamJSONFile::load(filename.c_str(), param))

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
  param.setValue("test:1:int_list", std::vector<int> {});
  param.setValue("test:1:double_list", std::vector<double> {});
  param.setValue("test:1:string_list", std::vector<std::string> {});
  param.setValue("test:1:file_output", std::string {}, "some description", {"output file"});
  param.setValue("test:1:is_executable_v1", std::string {}, "test is executable tag, giving a string", {"is_executable", "input file"});
  param.setValue("test:1:is_executable_v2", std::string {}, "test is executable tag, giving a type: File", {"is_executable", "input file"});


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
         "  \"string_list\": [\"SeqAn\", \"rocks\"],\n"
         "  \"file_output\": \"/some/made/up/path\",\n"
         "  \"is_executable_v1\": \"/some/made/up/path\",\n"
         "  \"is_executable_v2\": {\n"
         "        \"class\": \"File\",\n"
         "        \"path\": \"/some/made/up/path\"\n"
         "  }\n"
         "}\n";
  ofs.close();
  ParamJSONFile::load(filename.c_str(), param);

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
