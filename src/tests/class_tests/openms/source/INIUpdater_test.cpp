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
#include <OpenMS/APPLICATIONS/INIUpdater.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(INIUpdater, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

INIUpdater* ptr = nullptr;
INIUpdater* null_ptr = nullptr;
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
  TEST_EQUAL(m[Internal::ToolDescriptionInternal("FeatureFinder",ListUtils::create<String>("centroided"))]
             == Internal::ToolDescriptionInternal("FeatureFinderCentroided",ListUtils::create<String>("")), true)

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



