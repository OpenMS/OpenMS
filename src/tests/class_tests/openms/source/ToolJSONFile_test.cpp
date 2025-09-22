// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ToolJSONFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ToolJSONFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ToolJSONFile* ptr = nullptr;
ToolJSONFile* null_ptr = nullptr;

START_SECTION(ToolJSONFile())
{
  // Static class, no constructor test needed
  TEST_TRUE(true)
}
END_SECTION

START_SECTION(~ToolJSONFile())
{
  // Static class, no destructor test needed  
  TEST_TRUE(true)
}
END_SECTION

START_SECTION((static String getDefaultConfigPath()))
{
  String path = ToolJSONFile::getDefaultConfigPath();
  TEST_NOT_EQUAL(path, String())
  TEST_TRUE(path.hasSuffix("tools.json"))
}
END_SECTION

START_SECTION((static bool load(const String& filename, std::vector<Internal::ToolDescription>& tools, std::map<String, String>& categories)))
{
  // Create a temporary JSON file for testing
  String tmp_file = "/tmp/test_tools.json";
  
  // Write a simple test JSON
  std::ofstream ofs(tmp_file);
  ofs << R"({
    "categories": {
      "cat_test": "Test Category"
    },
    "tools": {
      "TestTool": {
        "category": "cat_test"
      }
    }
  })";
  ofs.close();
  
  std::vector<Internal::ToolDescription> tools;
  std::map<String, String> categories;
  
  bool result = ToolJSONFile::load(tmp_file, tools, categories);
  TEST_TRUE(result)
  TEST_EQUAL(tools.size(), 1)
  TEST_EQUAL(tools[0].name, "TestTool")
  TEST_EQUAL(tools[0].category, "Test Category")
  TEST_EQUAL(categories.size(), 1)
  TEST_EQUAL(categories["cat_test"], "Test Category")
  
  // Clean up
  std::remove(tmp_file.c_str());
}
END_SECTION

START_SECTION(Exception testing)
{
  std::vector<Internal::ToolDescription> tools;
  std::map<String, String> categories;
  
  // Test with non-existent file
  TEST_EXCEPTION(Exception::FileNotFound, ToolJSONFile::load("/non/existent/file.json", tools, categories))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST