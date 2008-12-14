// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SYSTEM/ExternalAllocatorUnique.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExternalAllocatorUnique, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExternalAllocatorUnique* ptr = 0;
String filename;
NEW_TMP_FILE(filename);

START_SECTION(ExternalAllocatorUnique())
{
	ptr = new ExternalAllocatorUnique(filename, 999);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ExternalAllocatorUnique())
{
	delete ptr;
}
END_SECTION

START_SECTION((ExternalAllocatorUnique(const String &filename, const Int64 &filesize)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	NOT_TESTABLE
}
END_SECTION

START_SECTION((ExternalAllocatorUnique(const ExternalAllocatorUnique &rhs)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	ExternalAllocatorUnique eau2(eau);
	NOT_TESTABLE
}
END_SECTION

START_SECTION((const String& getFilename() const))
{
	String filename2;
	NEW_TMP_FILE(filename2);
 	ExternalAllocatorUnique eau(filename2, 10000);
	TEST_EQUAL(eau.getFilename(), filename2)
}
END_SECTION

START_SECTION((const Int64& getFilesize() const))
{
  ExternalAllocatorUnique eau(filename, 10000);
	TEST_EQUAL(eau.getFilesize(), 10000)
}
END_SECTION

START_SECTION((void advanceFilesize(const Int64 &x)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	eau.advanceFilesize(33);
	eau.advanceFilesize(11);
	TEST_EQUAL(eau.getFilesize(), 10000+44);
}
END_SECTION

#ifdef OPENMS_WINDOWSPLATFORM
START_SECTION([EXTRA](const HANDLE& getMmapHandle() const))
#else
START_SECTION([EXTRA](const int& getMmapHandle() const))
#endif
{
	ExternalAllocatorUnique eau(filename, 10000);
  #ifdef OPENMS_WINDOWSPLATFORM
	eau.getMmapHandle();
	#else
	eau.getMmapHandle();
	#endif
	//hard to see if the handle is correct...
	NOT_TESTABLE
}
END_SECTION

START_SECTION((const Int64& getNextfree() const))
{
	ExternalAllocatorUnique eau(filename, 10000);
	TEST_EQUAL(eau.getNextfree(), 0);
}
END_SECTION

START_SECTION((void advanceNextfree(const Int64 &x)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	eau.advanceNextfree(33);
	eau.advanceNextfree(11);
	TEST_EQUAL(eau.getNextfree(), 44);
}
END_SECTION

START_SECTION((const Int64& getTotalmappingsize() const))
{
	ExternalAllocatorUnique eau(filename, 10000);
	TEST_EQUAL(eau.getTotalmappingsize(), 0);
}
END_SECTION

START_SECTION((void setTotalmappingsize(const Int64 &x)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	eau.setTotalmappingsize(33);
	TEST_EQUAL(eau.getTotalmappingsize(), 33);
}
END_SECTION

START_SECTION((bool hasFreeSwap(const Int64& bytes_needed)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	eau.advanceNextfree(33);
	TEST_EQUAL(eau.hasFreeSwap(9900), true);
	TEST_EQUAL(eau.hasFreeSwap(9990), false);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



