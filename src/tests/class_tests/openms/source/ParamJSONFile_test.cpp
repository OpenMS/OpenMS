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
#ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshadow"
#endif

using namespace OpenMS;

START_TEST(ParamJSONFile, "$Id")

START_SECTION((bool ParamJSONFile::load(const std::string& filename, Param& param)))
{
  std::string filename;
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

START_SECTION([EXTRA] bool ParamJSONFile::load() reads every JSON shape of a file parameter)
{
  // A CWL runner serializes a 'File[]' input as an array of 'File' objects. That shape used to
  // throw a json::type_error, which made every TOPP tool with an input file list unusable through
  // its generated CWL description as soon as the list was set (issue #10121).
  std::string filename;
  NEW_TMP_FILE(filename)

  // a fresh parameter tree per shape, so a value can never be left over from a previous load
  auto makeParam = []() {
    Param param;
    param.setValue("test:1:in", std::vector<std::string> {}, "input file list", {"input file"});
    param.setValue("test:1:out", std::vector<std::string> {}, "output file list", {"output file"});
    param.setValue("test:1:single", std::string {}, "single input file", {"input file"});
    param.setValue("test:1:plain_list", std::vector<std::string> {}, "untagged string list");
    return param;
  };
  auto loadJSON = [&filename](const std::string& content, Param& param) {
    std::ofstream ofs(filename.c_str(), std::ios::out);
    ofs << content;
    ofs.close();
    ParamJSONFile::load(filename.c_str(), param);
  };

  // (a) array of CWL 'File' objects -- what a CWL runner writes for 'type: File[]'
  Param param_a = makeParam();
  loadJSON(R"({"in": [{"class": "File", "path": "a.mzML"}, {"class": "File", "path": "b.mzML"}]})", param_a);
  std::vector<std::string> files_a = param_a.getValue("test:1:in").toStringVector();
  TEST_EQUAL(files_a.size(), 2);
  TEST_STRING_EQUAL(files_a[0], "a.mzML");
  TEST_STRING_EQUAL(files_a[1], "b.mzML");

  // (b) array of plain path strings
  Param param_b = makeParam();
  loadJSON(R"({"in": ["a.mzML", "b.mzML"]})", param_b);
  std::vector<std::string> files_b = param_b.getValue("test:1:in").toStringVector();
  TEST_EQUAL(files_b.size(), 2);
  TEST_STRING_EQUAL(files_b[0], "a.mzML");
  TEST_STRING_EQUAL(files_b[1], "b.mzML");

  // (c) object with a 'path' array -- what ParamJSONFile::store writes, with and without 'class'
  Param param_c = makeParam();
  loadJSON(R"({"in": {"class": "File", "path": ["a.mzML", "b.mzML"]}})", param_c);
  std::vector<std::string> files_c = param_c.getValue("test:1:in").toStringVector();
  TEST_EQUAL(files_c.size(), 2);
  TEST_STRING_EQUAL(files_c[0], "a.mzML");
  TEST_STRING_EQUAL(files_c[1], "b.mzML");

  Param param_c2 = makeParam();
  loadJSON(R"({"in": {"path": ["a.mzML", "b.mzML"]}})", param_c2);
  std::vector<std::string> files_c2 = param_c2.getValue("test:1:in").toStringVector();
  TEST_EQUAL(files_c2.size(), 2);
  TEST_STRING_EQUAL(files_c2[0], "a.mzML");
  TEST_STRING_EQUAL(files_c2[1], "b.mzML");

  // a single file for a list-valued parameter, in both notations
  Param param_single_as_list = makeParam();
  loadJSON(R"({"in": {"class": "File", "path": "a.mzML"}})", param_single_as_list);
  std::vector<std::string> files_single = param_single_as_list.getValue("test:1:in").toStringVector();
  TEST_EQUAL(files_single.size(), 1);
  TEST_STRING_EQUAL(files_single[0], "a.mzML");

  // 'Directory' objects carry their path the same way 'File' objects do
  Param param_dir = makeParam();
  loadJSON(R"({"in": [{"class": "Directory", "path": "/data/run1.d"}]})", param_dir);
  std::vector<std::string> dirs = param_dir.getValue("test:1:in").toStringVector();
  TEST_EQUAL(dirs.size(), 1);
  TEST_STRING_EQUAL(dirs[0], "/data/run1.d");

  // an empty list stays empty
  Param param_empty = makeParam();
  loadJSON(R"({"in": []})", param_empty);
  TEST_EQUAL(param_empty.getValue("test:1:in").toStringVector().size(), 0);

  // output file lists are CWL strings, but accept the same shapes
  Param param_out = makeParam();
  loadJSON(R"({"out": ["a.mzML", "b.mzML"], "in": [{"class": "File", "path": "c.mzML"}]})", param_out);
  std::vector<std::string> out_files = param_out.getValue("test:1:out").toStringVector();
  TEST_EQUAL(out_files.size(), 2);
  TEST_STRING_EQUAL(out_files[0], "a.mzML");
  TEST_STRING_EQUAL(out_files[1], "b.mzML");

  // a single input file, as an object and as a plain string
  Param param_scalar = makeParam();
  loadJSON(R"({"single": {"class": "File", "path": "a.mzML"}})", param_scalar);
  TEST_STRING_EQUAL(std::string(param_scalar.getValue("test:1:single")), "a.mzML");

  Param param_scalar_str = makeParam();
  loadJSON(R"({"single": "a.mzML"})", param_scalar_str);
  TEST_STRING_EQUAL(std::string(param_scalar_str.getValue("test:1:single")), "a.mzML");

  // an untagged string list is unaffected by all of this
  Param param_plain = makeParam();
  loadJSON(R"({"plain_list": ["SeqAn", "rocks"]})", param_plain);
  TEST_EQUAL(param_plain.getValue("test:1:plain_list").toStringVector().size(), 2);

  // shapes that carry no path at all are rejected, rather than read as an empty list
  Param param_bad = makeParam();
  TEST_EXCEPTION(Exception::ParseError, loadJSON(R"({"in": 42})", param_bad))
  TEST_EXCEPTION(Exception::ParseError, loadJSON(R"({"in": [{"class": "File"}]})", param_bad))
  TEST_EXCEPTION(Exception::ParseError, loadJSON(R"({"in": {"class": "File"}})", param_bad))
  TEST_EXCEPTION(Exception::ParseError, loadJSON(R"({"single": 42})", param_bad))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#ifdef __clang__
  #pragma clang diagnostic pop
#endif