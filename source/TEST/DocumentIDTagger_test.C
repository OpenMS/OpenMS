// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <cstdio>
#include <cstdlib>
///////////////////////////
#include <OpenMS/METADATA/DocumentIDTagger.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DocumentIDTagger, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DocumentIDTagger* ptr = 0;
START_SECTION(DocumentIDTagger())
{
	ptr = new DocumentIDTagger("someTOPPTool");
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DocumentIDTagger())
{
	delete ptr;
}
END_SECTION

START_SECTION((DocumentIDTagger(String toolname)))
{
	DocumentIDTagger tagme("SomeTOPPTool");
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DocumentIDTagger(const DocumentIDTagger &source)))
{
  DocumentIDTagger tagme("SomeTOPPTool");
  DocumentIDTagger tagme2(tagme);

  TEST_EQUAL(tagme==tagme2, true)
}
END_SECTION

START_SECTION((DocumentIDTagger& operator=(const DocumentIDTagger &source)))
{
  DocumentIDTagger tagme("SomeTOPPTool");
	DocumentIDTagger tagme2 = tagme;
	TEST_EQUAL(tagme==tagme2,true)
}
END_SECTION

START_SECTION((bool operator==(const DocumentIDTagger &source) const ))
{
  DocumentIDTagger tagme("SomeTOPPTool");
	DocumentIDTagger tagme2 = tagme;
	TEST_EQUAL(tagme==tagme2, true)
	DocumentIDTagger tagme3(tagme);
	TEST_EQUAL(tagme==tagme3, true)
}
END_SECTION

START_SECTION((bool operator!=(const DocumentIDTagger &source) const ))
{
  DocumentIDTagger tagme("SomeTOPPTool");
	DocumentIDTagger tagme2("SomeOtherTOPPTool");
	TEST_EQUAL(tagme!=tagme2, true)
}
END_SECTION

START_SECTION((String getPoolFile() const))
{
	NOT_TESTABLE; // tested below
}
END_SECTION
			
START_SECTION((void setPoolFile(const String& file)))
{
	String tmp_pool;
	NEW_TMP_FILE(tmp_pool);
	DocumentIDTagger tagme("SomeTOPPTool");
	//use custom pool file
	tagme.setPoolFile(tmp_pool);
	TEST_EQUAL(tagme.getPoolFile(), tmp_pool)
}
END_SECTION			

String tmp_pool;
NEW_TMP_FILE(tmp_pool);
ofstream outfile;
outfile.open (tmp_pool.c_str(), ofstream::out);
outfile << "ID1\nIDNew\nIDsecondtoLast\nIDLast\n";
outfile.close();

START_SECTION((bool tag(DocumentIdentifier &map) const ))
{
	DocumentIdentifier myD;
	myD.setIdentifier("");
	DocumentIDTagger tagme("SomeTOPPTool");
	//use custom pool file
	tagme.setPoolFile(tmp_pool);
	Int cnt(0);
	tagme.countFreeIDs(cnt);
	TEST_EQUAL(cnt, 4);
	tagme.tag(myD);
	TEST_EQUAL(myD.getIdentifier(), "ID1");
	tagme.tag(myD);
	TEST_EQUAL(myD.getIdentifier(), "IDNew");
	tagme.countFreeIDs(cnt);
	TEST_EQUAL(cnt, 2);
	// 2 left
	TEST_EQUAL(tagme.tag(myD), true);
	// 1 left
	TEST_EQUAL(tagme.tag(myD), true);
	//0 left, expect it to go wrong
	TEST_EXCEPTION(Exception::DepletedIDPool, tagme.tag(myD));
	// confirm 0 left
	TEST_EQUAL(tagme.countFreeIDs(cnt), true);
	TEST_EQUAL(cnt, 0);
}
END_SECTION

START_SECTION((bool countFreeIDs(Int &free) const ))
{
  NOT_TESTABLE //done above
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



