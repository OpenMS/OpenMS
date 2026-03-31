// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ToolDescriptionFile.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
///////////////////////////

#include <filesystem>
#include <vector>
#include <string>

using namespace OpenMS;
using namespace std;

START_TEST(ToolDescriptionFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ToolDescriptionFile* ptr = nullptr;
ToolDescriptionFile* null_ptr = nullptr;
START_SECTION(ToolDescriptionFile())
{
	ptr = new ToolDescriptionFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~ToolDescriptionFile())
{
	delete ptr;
}
END_SECTION

START_SECTION((void load(const String &filename, std::vector< Internal::ToolDescription > &tds)))
{
  ToolDescriptionFile f;
  std::vector< Internal::ToolDescription > tds;
  std::filesystem::path dir_path{std::string(ToolHandler::getExternalToolsPath())};
  std::vector<std::string> files;
  for (const auto& entry : std::filesystem::directory_iterator(dir_path))
  {
    if (entry.path().extension() == ".ttd")
    {
      files.push_back(entry.path().string());
    }
  }
  for (const auto& file : files)
  {
    f.load(file, tds);
    //std::cerr << "load: " << file << "\n";
    TEST_EQUAL(!tds.empty(), true)
  }
  
}
END_SECTION

START_SECTION((void store(const String &filename, const std::vector< Internal::ToolDescription > &tds) const ))
{
  ToolDescriptionFile f;
  std::vector< Internal::ToolDescription > tds;
  TEST_EXCEPTION( Exception::NotImplemented, f.store("bla", tds))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



