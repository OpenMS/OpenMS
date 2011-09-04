// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/ToolHandler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ToolHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ToolHandler* ptr = 0;
ToolHandler* null_ptr = 0;
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
  TEST_EQUAL(list.has("FeatureFinderMRM"), true)
  TEST_EQUAL(list.has("GenericWrapper"), false)
  TEST_EQUAL(list.size() > 30, true)  // assume we have over 30 tools in there
  list = ToolHandler::getTOPPToolList(true);
  TEST_EQUAL(list.has("FeatureFinderMRM"), true)
  TEST_EQUAL(list.has("GenericWrapper"), true)
  TEST_EQUAL(list.size() > 30, true)  // assume we have over 30 tools in there
}
END_SECTION

START_SECTION((static ToolListType getUtilList()))
{
  ToolListType list = ToolHandler::getUtilList();
  TEST_EQUAL(list.has("ImageCreator"), true)
  TEST_EQUAL(list.has("FFEval"), true)
  TEST_EQUAL(list.size() > 10, true)  // assume we have over 30 tools in there
}
END_SECTION

START_SECTION((static StringList getTypes(const String &toolname)))
{
  StringList  results = StringList::create("4plex,8plex");
  TEST_EQUAL(ToolHandler::getTypes("ITRAQAnalyzer"), results);

  TEST_EQUAL(ToolHandler::getTypes("IDMapper"), StringList());
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



