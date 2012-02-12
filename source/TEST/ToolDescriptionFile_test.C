// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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

ToolDescriptionFile* ptr = 0;
ToolDescriptionFile* null_ptr = 0;
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
    TEST_EQUAL(tds.size()>=1, true)
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



