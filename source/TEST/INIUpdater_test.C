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
#include <OpenMS/APPLICATIONS/INIUpdater.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(INIUpdater, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

INIUpdater* ptr = 0;
INIUpdater* null_ptr = 0;
START_SECTION(INIUpdater())
{
	ptr = new INIUpdater();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~INIUpdater())
{
	delete ptr;
}
END_SECTION

START_SECTION((StringList getToolNamesFromINI(const Param &ini) const))
{
  Param p;
  INIUpdater i;
  StringList names = i.getToolNamesFromINI(p);

  TEST_EQUAL(names.size(), 0)

  p.setValue("FeatureFinder:version","1.9");
  p.setValue("SomeTool:version","whatever");
  names = i.getToolNamesFromINI(p);

  TEST_EQUAL(names.size(), 2)

  p.setValue("BrokenTool:version2","1.9");
  names = i.getToolNamesFromINI(p);

  TEST_EQUAL(names.size(), 2)

}
END_SECTION

START_SECTION((const ToolMapping& getNameMapping()))
{
  INIUpdater i;
  ToolMapping m = i.getNameMapping();

  TEST_NOT_EQUAL(m.size(), 0)
  TEST_EQUAL(m[Internal::ToolDescriptionInternal("FeatureFinder",StringList::create("centroided"))]
             == Internal::ToolDescriptionInternal("FeatureFinderCentroided",StringList::create("")), true)

}
END_SECTION

START_SECTION((bool getNewToolName(const String &old_name, const String &tools_type, String &new_name)))
{
  INIUpdater i;
  String new_name;
  i.getNewToolName("FeatureFinder", "centroided", new_name);
  TEST_EQUAL(new_name, "FeatureFinderCentroided");

  i.getNewToolName("PeakPicker", "wavelet", new_name);
  TEST_EQUAL(new_name, "PeakPickerWavelet");

  i.getNewToolName("ITRAQAnalyzer", "4plex", new_name);
  TEST_EQUAL(new_name, "ITRAQAnalyzer");

  i.getNewToolName("FileInfo", "", new_name);
  TEST_EQUAL(new_name, "FileInfo");

  i.getNewToolName("FileInfo", "bogus type", new_name); // type will be ignored - ok
  TEST_EQUAL(new_name, "FileInfo");

  TEST_EQUAL(i.getNewToolName("UNKNOWNTOOL", "bogus type", new_name), false);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



