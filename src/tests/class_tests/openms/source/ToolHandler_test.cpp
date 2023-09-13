// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

START_SECTION((static ToolListType getTOPPToolList(const bool includeGenericWrapper=false)))
{
  ToolListType list = ToolHandler::getTOPPToolList();
  TEST_EQUAL(list.find("FeatureFinderMRM") != list.end(), true)
  TEST_EQUAL(list.find("GenericWrapper") != list.end(), false)
  TEST_EQUAL(list.size() > 30, true)  // assume we have over 30 tools in there
  list = ToolHandler::getTOPPToolList(true);
  TEST_EQUAL(list.find("FeatureFinderMRM") != list.end(), true)
  TEST_EQUAL(list.find("GenericWrapper") != list.end(), true)
  TEST_EQUAL(list.size() > 30, true)  // assume we have over 30 tools in there
}
END_SECTION

START_SECTION((static StringList getTypes(const String &toolname)))
{
  TEST_EQUAL(ToolHandler::getTypes("IsobaricAnalyzer").empty(), true);
  TEST_EQUAL(ToolHandler::getTypes("IDMapper").empty(), true);
}
END_SECTION

START_SECTION((static String getExternalToolsPath()))
{
  TEST_NOT_EQUAL(ToolHandler::getExternalToolsPath(), String())
}
END_SECTION

START_SECTION((static String getInternalToolsPath()))
{
  TEST_NOT_EQUAL(ToolHandler::getExternalToolsPath(), String())
}
END_SECTION

START_SECTION((static String getCategory(const String &toolname)))
{
  TEST_EQUAL(ToolHandler::getCategory("IDFilter"), "ID Processing")
  TEST_EQUAL(ToolHandler::getCategory("DOESNOTEXIST"), "")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



