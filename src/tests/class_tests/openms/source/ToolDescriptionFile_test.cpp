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

#include <QStringList>
#include <QDir>

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
  QStringList list;
  QDir dir( ToolHandler::getExternalToolsPath().toQString(), "*.ttd");
  QStringList files = dir.entryList();
  for (int i=0;i<files.size();++i)
  {
    files[i] = dir.absolutePath()+QDir::separator()+files[i];
    f.load(files[i], tds);
    //std::cerr << "load: " << String(files[i]) << "\n";
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



