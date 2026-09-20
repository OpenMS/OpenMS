// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/ToolHandler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ToolHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ToolHandler* ptr = nullptr;
ToolHandler* null_ptr = nullptr;
START_SECTION(ToolHandler())
{
	ptr = new ToolHandler();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ToolHandler())
{
	delete ptr;
}
END_SECTION

START_SECTION((static ToolListType getTOPPToolList()))
{
  ToolListType list = ToolHandler::getTOPPToolList();
  TEST_TRUE(list.find("DecoyDatabase") != list.end())
  TEST_TRUE(list.size() > 30)  // assume we have over 30 tools in there
  // the registry comes from the registry files of the share directory; an empty list means they
  // were not found, which would let every check below pass vacuously
  TEST_TRUE(!list.empty())
  // the key and the description's own name are the same thing: the name the registry entry declares
  TEST_EQUAL(list.find("DecoyDatabase")->second.name, "DecoyDatabase")
  TEST_EQUAL(list.find("QCShrinker")->second.name, "QCShrinker")
  // a tool of a build option is listed only by a build that has the option, because
  // ToolHandler reads that option's registry file only then
#ifdef WITH_GUI
  TEST_TRUE(list.find("ImageCreator") != list.end())
#else
  TEST_TRUE(list.find("ImageCreator") == list.end())
#endif
#ifdef WITH_WNETALIGN
  TEST_TRUE(list.find("FeatureLinkerWNet") != list.end())
#else
  TEST_TRUE(list.find("FeatureLinkerWNet") == list.end())
#endif
}
END_SECTION

START_SECTION((static const ToolListType& getTOPPToolListRef()))
{
  const ToolListType& list = ToolHandler::getTOPPToolListRef();
  TEST_EQUAL(list.size(), ToolHandler::getTOPPToolList().size())
  TEST_TRUE(list.find("DecoyDatabase") != list.end())
  // served from the process-wide cache, so the reference is the same one every time
  TEST_EQUAL(&list, &ToolHandler::getTOPPToolListRef())
}
END_SECTION

START_SECTION((static StringList getTypes(const std::string &toolname)))
{
  TEST_EQUAL(ToolHandler::getTypes("IsobaricAnalyzer").empty(), true);
  TEST_EQUAL(ToolHandler::getTypes("IDMapper").empty(), true);
  // An unknown tool has no types rather than being an error: a tool that is not in this
  // installation's registry still has to be able to write its own CTD/CWL description
  // (TOPPBase::handleWriteCommands_ asks for the types to know how many files to write).
  TEST_EQUAL(ToolHandler::getTypes("DOESNOTEXIST").empty(), true);
}
END_SECTION

START_SECTION((static std::string getToolRegistryPath()))
{
  TEST_NOT_EQUAL(ToolHandler::getToolRegistryPath(), std::string())
}
END_SECTION

START_SECTION((static std::string getCategory(const std::string &toolname)))
{
  TEST_EQUAL(ToolHandler::getCategory("IDFilter"), "File Filtering, Extraction and Merging")
  TEST_EQUAL(ToolHandler::getCategory("DOESNOTEXIST"), "")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
