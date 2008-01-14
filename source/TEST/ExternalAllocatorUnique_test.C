// -*- Mode: C++; tab-width: 2; -*-
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

CHECK(ExternalAllocatorUnique())
{
	ptr = new ExternalAllocatorUnique(filename, 999);
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~ExternalAllocatorUnique())
{
	delete ptr;
}
RESULT

CHECK((ExternalAllocatorUnique(const String &filename, const Offset64Int &filesize)))
{
	ExternalAllocatorUnique eau(filename, 10000);
}
RESULT

CHECK((ExternalAllocatorUnique(const ExternalAllocatorUnique &rhs)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	ExternalAllocatorUnique eau2(eau);
}
RESULT

CHECK((const String& getFilename() const))
{
	String filename2;
	NEW_TMP_FILE(filename2);
 	ExternalAllocatorUnique eau(filename2, 10000);
	TEST_EQUAL(eau.getFilename(), filename2)
}
RESULT

CHECK((const Offset64Int& getFilesize() const))
{
  ExternalAllocatorUnique eau(filename, 10000);
	TEST_EQUAL(eau.getFilesize(), 10000)
}
RESULT


#ifdef OPENMS_WINDOWSPLATFORM
CHECK((const HANDLE& getMmapHandle() const))
#else
CHECK((const int& getMmapHandle() const))
#endif
{
	ExternalAllocatorUnique eau(filename, 10000);
  #ifdef OPENMS_WINDOWSPLATFORM
	HANDLE h = eau.getMmapHandle();
	#else
	int h = eau.getMmapHandle();
	#endif
	//hard to see if the handle is correct...
}
RESULT

CHECK((const Offset64Int& getNextfree() const))
{
	ExternalAllocatorUnique eau(filename, 10000);
	TEST_EQUAL(eau.getNextfree(), 0);
}
RESULT

CHECK((void advanceNextfree(const Offset64Int &x)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	eau.advanceNextfree(33);
	eau.advanceNextfree(11);
	TEST_EQUAL(eau.getNextfree(), 44);
}
RESULT

CHECK((const Offset64Int& getTotalmappingsize() const))
{
	ExternalAllocatorUnique eau(filename, 10000);
	TEST_EQUAL(eau.getTotalmappingsize(), 0);
}
RESULT

CHECK((void setTotalmappingsize(const Offset64Int &x)))
{
	ExternalAllocatorUnique eau(filename, 10000);
	eau.setTotalmappingsize(33);
	TEST_EQUAL(eau.getTotalmappingsize(), 33);
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



